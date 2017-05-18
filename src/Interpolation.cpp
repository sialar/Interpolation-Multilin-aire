#include "../include/Interpolation.hpp"

Interpolation::~Interpolation()
{}
Interpolation::Interpolation(int d, int nIter, int method)
{
    m_method = method;
    m_d = d;
    m_maxIteration = nIter;

    m_interpolationPoints.resize(m_d);
    m_middles.resize(m_d);

    m_lejaSequence = Utils::loadLejaSequenceFromFile(m_maxIteration);
    m_middlePoints = Utils::createSequenceByDichotomy(m_maxIteration);

    for (int i=0; i<m_d; i++)
        if (m_method == 1)
            m_middles[i] = new BinaryTree(int(log(m_maxIteration)/log(2)));
}

/************************* Data Points ****************************************/
MultiVariatePoint<double> Interpolation::getPoint(MultiVariatePoint<int> nu)
{
    MultiVariatePoint<double> point(m_d,0.0);
    for (int i=0; i<m_d; i++)
        if (m_method == 0)
            point(i) = m_lejaSequence[nu(i)];
        else
            point(i) = m_middlePoints[nu(i)];
    return point;
}
void Interpolation::addInterpolationPoint(MultiVariatePoint<double> p)
{
    m_interpolationNodes.push_back(p);

    bool found;
    for (int i=0; i<m_d; i++)
    {
        found = false;
        for (double x : m_interpolationPoints[i])
            if (x == p(i)) found = true;
        if (!found) m_interpolationPoints[i].push_back(p(i));
    }
}

/******************************************************************************/


/************************* AI algo ********************************************/
int Interpolation::buildPathWithAIAlgo(auto start_time, double threshold, bool debug)
{
    m_curentNeighbours.clear();
    m_path.clear();

    MultiVariatePoint<int>* firstPoint = new MultiVariatePoint<int>(m_d,0);
    m_path.push_back(firstPoint);
    addInterpolationPoint(getPoint(*firstPoint));
    m_alphaMap.insert(pair<MultiVariatePoint<int>*,double>(firstPoint,Utils::gNd(getPoint(*firstPoint))));

    double val, max;
    MultiVariatePoint<int>* argmax = new MultiVariatePoint<int>(m_d,0);
    m_curentNeighbours.push_back(argmax);
    int iteration = 1;

    while (!m_curentNeighbours.empty() && iteration < m_maxIteration)
    {
        max = 0;
        if (m_method)
            updateNextPoints(argmax);
        else
            updateCurentNeighbours(argmax);

        if (debug)
        {
            Utils::separateur();
            displayPath();
            displayCurentNeighbours();
            displayAlphaTab();
        }


        for (MultiVariatePoint<int>* nu : m_curentNeighbours)
        {
            map<MultiVariatePoint<int>*,double>::iterator it = m_alphaMap.find(nu);
            if (it != m_alphaMap.end())
                val = m_alphaMap[nu];
            else
                val =  computeLastAlphaNu(nu);

            nu->setAlpha(val);
            if (debug) cout << *nu << " " << val << " | ";

            if (abs(val) >= abs(max))
            {
                argmax = nu;
                max = val;
            }
        }

        double e = 1e-10;
        bool allAlphaLowerThanE = true;
        for (MultiVariatePoint<int>* nu : m_curentNeighbours)
            if (nu->getAlpha()>e) allAlphaLowerThanE = false;

        if (allAlphaLowerThanE)
        {
            int maxAge = -1;
            for (MultiVariatePoint<int>* nu : m_curentNeighbours)
            {
                if (nu->getWaitingTime()>maxAge)
                {
                    maxAge = nu->getWaitingTime();
                    argmax = nu;
                    max = nu->getAlpha();
                }
            }
        }

        m_path.push_back(argmax);
        m_alphaMap.insert(pair<MultiVariatePoint<int>*,double>(argmax,max));

        addInterpolationPoint(getPoint(*argmax));
        iteration++;
        if (debug) cout << endl;

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI

        if ((m_maxIteration>10) && iteration%(m_maxIteration/10)==0)
        {
            auto end_time = chrono::steady_clock::now();
            std::chrono::duration<double> run_time = end_time - start_time;

            double error = tryWithCurentPath();
            cout << "   - Interpolation error after " << iteration << " iterations: " << error;
            cout << " | Elapsed time : "  << run_time.count() << endl;
            if (error < threshold)
            {
                cout << endl << "   - AI Algo stop after " << iteration << " iterations";
                cout << " | Elapsed time : "  << run_time.count() << endl;
                return iteration;
            }
        }

    }
    cout << "Number of iterations : " << iteration << endl;
    return iteration;
}
void Interpolation::updateCurentNeighbours(MultiVariatePoint<int>* nu)
{
  m_curentNeighbours.remove(nu);
  for (MultiVariatePoint<int>* mu : m_curentNeighbours)
      mu->incrWaitingTime();

  for (int k=0; k<nu->getD(); k++)
  {
      MultiVariatePoint<int>* nu1 = new MultiVariatePoint<int>(*nu);
      (*nu1)(k)++;
      if (isCorrectNeighbourToCurentPath(nu1))
          m_curentNeighbours.push_back(nu1);
  }
}
void Interpolation::updateNextPoints(MultiVariatePoint<int>* nu)
{
    m_curentNeighbours.remove(nu);
    for (MultiVariatePoint<int>* mu : m_curentNeighbours)
        mu->incrWaitingTime();

    for (int i=0; i<m_d; i++)
    {
        double inf, sup;
        Node* node = m_middles[i]->searchNode(getPoint(*nu)(i),&inf,&sup,true);
        if (node->left())
        {
            MultiVariatePoint<int>* mu1 = new MultiVariatePoint<int>(*nu);
            inf = node->left()->key();
            (*mu1)(i) = BinaryTree::getIndice(inf);
            if (!indiceInPath(*mu1) && !indiceInNeighborhood(*mu1)) m_curentNeighbours.push_back(mu1);
        }
        if (node->right())
        {
            MultiVariatePoint<int>* mu2 = new MultiVariatePoint<int>(*nu);
            sup = node->right()->key();
            (*mu2)(i) = BinaryTree::getIndice(sup);
            if (!indiceInPath(*mu2) && !indiceInNeighborhood(*mu2) ) m_curentNeighbours.push_back(mu2);
        }
    }
}
bool Interpolation::isCorrectNeighbourToCurentPath(MultiVariatePoint<int>* nu)
{
    bool isNeighbour[m_d];
    for (int i=0; i<m_d; i++)
        isNeighbour[i] = true;

    MultiVariatePoint<int> index(*nu);
    for (int i=0; i<m_d; i++)
        if ((*nu)(i))
        {
            index(i)--;
            isNeighbour[i] = indiceInPath(index);
            index(i)++;
        }
    bool res = true;
    for (int i=0; i<m_d; i++)
        res = res && isNeighbour[i];
    return res;
}
bool Interpolation::indiceInNeighborhood(MultiVariatePoint<int> index)
{
    for (MultiVariatePoint<int>* nu : m_curentNeighbours)
        if (index==*nu)
            return true;
    return false;
}
bool Interpolation::indiceInPath(MultiVariatePoint<int> index)
{
    for (MultiVariatePoint<int>* nu : m_path)
        if (index==*nu)
            return true;
    return false;
}
double Interpolation::computeLastAlphaNu(MultiVariatePoint<int>* nu)
{
    double basisFuncProd;
    double res = Utils::gNd(getPoint(*nu));
    for (MultiVariatePoint<int>* l : m_path)
    {
        basisFuncProd = 1;
        for (int p=0; p<m_d; p++)
        {
            if (m_method == 0)
                basisFuncProd *= lagrangeBasisFunction_1D((*l)(p),getPoint(*nu)(p),p);
            else if (m_method == 1)
                basisFuncProd *= piecewiseFunction_1D((*l)(p),getPoint(*nu)(p),p); // TODO: getPoint() for method 1
        }
        res -= m_alphaMap[l] * basisFuncProd;
    }
    return res;
}
/******************************************************************************/


/*********************** Test functions ***************************************/
void Interpolation::testPathBuilt(double threshold, bool debug)
{
  auto start_time = chrono::steady_clock::now();
  int nbIterations = buildPathWithAIAlgo(start_time, threshold, debug);
  auto end_time = chrono::steady_clock::now();
  std::chrono::duration<double> run_time = end_time - start_time;
  cout << endl << " - Time required to compute the path with " << nbIterations <<
  " iterations: " << run_time.count() << "s" << endl;
}
double Interpolation::tryWithCurentPath()
{
  vector<double> realValues, estimate;
  for (MultiVariatePoint<double> p : m_testPoints)
  {
    realValues.push_back(Utils::gNd(p));
    estimate.push_back(interpolation_ND(p));
  }
  return Utils::interpolationError(realValues,estimate);
}
/******************************************************************************/


/************************* Display functions **********************************/
void Interpolation::displayPath()
{
    // Format: (nu:points[nu]) --> () ...
    cout << "Chemin =";
    int n = m_path.size();
    for (int i=0; i<n; i++)
    {
        if (i>0) cout << "\t";
        cout << " " << i << " :";
        cout << " [" << *m_path[i] << ":" << getPoint(*m_path[i]) << ":" << m_path[i] << "]";
        cout << " --> alpha(" << *m_path[i] << ") = " << m_alphaMap[m_path[i]] << endl;
    }
}
void Interpolation::displayAlphaTab()
{
    map<MultiVariatePoint<int>*,double>::iterator it;
    cout << "Values of alpha (" << m_alphaMap.size() << ") :" << endl;
    for (it=m_alphaMap.begin(); it!=m_alphaMap.end(); it++)
        cout << "(" << *(get<0>(*it)) << ":" << getPoint(*(get<0>(*it))) << ":" << get<0>(*it) << ") : " << get<1>(*it) <<  endl;
}
void Interpolation::displayCurentNeighbours()
{
    cout << "Curent neighbours (" << m_curentNeighbours.size() << ") = ";
    for (MultiVariatePoint<int>* nu : m_curentNeighbours)
        cout << "(" << (*nu) << ":" << getPoint(*nu) << ":" << nu << ") [" << nu << " " << nu->getWaitingTime() << "] | ";
    cout << endl;
}
void Interpolation::displayInterpolationPointsInEachDirection()
{
    vector<double>::iterator it;
    for (int i=0; i<m_d; ++i)
    {
        cout << " - " << m_interpolationPoints[i].size() << " points in direction " << i << " : { ";
        for (it=m_interpolationPoints[i].begin(); it!=m_interpolationPoints[i].end(); it++)
            cout << *it << " ";
        cout << "}" << endl;
    }
}
void Interpolation::displayInterpolationMultiVariatePoints()
{
    cout << " - Interpolation nodes: { ";
    for (MultiVariatePoint<double> x : m_interpolationNodes)
        cout << x << " ";
    cout << "}" << endl;
}
void Interpolation::savePathInFile()
{
  ofstream file("python/path.txt", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==2 || m_d==3)
    {
      cout << " - The path is saved in python/path.txt" << endl;
      file << m_d << endl;
      file << m_interpolationPoints[0].size() << " " <<  m_interpolationPoints[1].size() << " ";
      if (m_d == 2) file << "1" << endl;
      else file << m_interpolationPoints[2].size() << endl;
      file << m_path.size() << endl;
      for (MultiVariatePoint<int>* nu : m_path)
      {
        file << (*nu)(0) << " " << (*nu)(1) << " ";
        if (m_d == 2) file << "0" << endl;
        else file << (*nu)(2) << endl;
      }
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}
/******************************************************************************/


/*********************** Interpolation ****************************************/
void Interpolation::storeInterpolationFunctions()
{
    ofstream file("python/plot_function.txt", ios::out | ios::trunc);
    if(file)
    {
      if (m_d==1)
      {
          vector<double> x;
          int nbPoints = int(m_testPoints.size());
          for (int i=0; i<nbPoints; i++)
              x.push_back(m_testPoints[i](0));
          sort(x.begin(), x.end());
          file << m_path.size() << " " << nbPoints << endl;
          for (int j=0; j<nbPoints; j++)
          {
              file << x[j];
              for (MultiVariatePoint<int>* nu : m_path)
              {
                  if (m_method == 0)
                      file << " " << lagrangeBasisFunction_1D((*nu)(0),x[j],0);
                  else if (m_method == 1)
                      file << " " << piecewiseFunction_1D((*nu)(0),x[j],0);
              }
              file << " " << interpolation_ND(MultiVariatePoint<double>::toMonoVariatePoint(x[j]));
              file << " " <<  Utils::gNd(MultiVariatePoint<double>::toMonoVariatePoint(x[j]));
              file << endl;
          }
      }
      file.close();
    }
    else
    cerr << "Error while opening the file!" << endl;
}

double Interpolation::piecewiseFunction_1D(int k, double t, int axis)
{
    if (k==0) return 1;
    else if (k==1) return (t>0) ? 0 : -t;
    else if (k==2) return (t<0) ? 0 : t;
    else
    {
        double sup, inf, tk = m_middlePoints[k];
        m_middles[axis]->searchNode(tk,&sup,&inf, false);
        if (t <= tk && t >= inf) return ((t-inf)/(tk-inf));
        else if (t >= tk && t <= sup) return ((t-sup)/(tk-sup));
        else return 0;
    }
}
double Interpolation::lagrangeBasisFunction_1D(int k, double t, int axis)
{
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        prod *= (t-m_lejaSequence[i]) / (m_lejaSequence[k]-m_lejaSequence[i]) ;
    return prod;
}
double Interpolation::interpolation_ND(MultiVariatePoint<double> x)
{
    double l_prod, sum = 0;
    for (MultiVariatePoint<int>* nu : m_path)
    {
        l_prod = 1;
        for (int i=0; i<m_d; i++)
        {
            if (m_method == 0)
                l_prod *= lagrangeBasisFunction_1D((*nu)(i),x(i),i);
            else if (m_method == 1)
                l_prod *= piecewiseFunction_1D((*nu)(i),x(i),i);
        }
        sum += l_prod * m_alphaMap[nu];
    }
    return sum;
}
/******************************************************************************/
