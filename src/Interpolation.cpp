#include "../include/Interpolation.hpp"

Interpolation::~Interpolation()
{
    //for (MultiVariatePointPtr<int> p : m_path)
    //    delete p;

    //map<MultiVariatePointPtr<int>,double>::iterator it;
    //for (it = m_alphaMap.begin(); it != m_alphaMap.end(); it++)
        //delete get<0>(*it);

}

Interpolation::Interpolation(int d, int nIter, int method)
{
    m_method = method;
    m_d = d;
    m_maxIteration = nIter;

    m_interpolationPoints.resize(m_d);

    m_lejaSequence = Utils::loadLejaSequenceFromFile(m_maxIteration);
    m_middlePoints = Utils::createSequenceByDichotomy(m_maxIteration);
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
void Interpolation::addInterpolationPoint(MultiVariatePoint<double>p)
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
void Interpolation::computeBoundaries(double t, double* inf, double* sup, int axis)
{
    vector<double> higherPoints, lowerPoints;
    for (double x : m_interpolationPoints[axis])
    {
        if (x>t) higherPoints.push_back(x);
        else if (x<t) lowerPoints.push_back(x);
        else break;
    }
    if (higherPoints.size())
        *inf = *min_element(higherPoints.begin(),higherPoints.end());
    else *inf = t;
    if (lowerPoints.size())
        *sup = *max_element(lowerPoints.begin(),lowerPoints.end());
    else *sup = t;
}

/******************************************************************************/


/************************* AI algo ********************************************/
bool alphaLess(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    return abs(nu->getAlpha()) < abs(mu->getAlpha());
}
bool ageLess(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    return abs(nu->getWaitingTime()) < abs(mu->getWaitingTime());
}
int Interpolation::buildPathWithAIAlgo(auto start_time, double threshold, bool debug)
{
    m_curentNeighbours.clear();
    m_path.clear();

    int iteration = 0;
    double val, max, e = 1e-10;
    bool allAlphaLowerThanE = true;
    MultiVariatePointPtr<int> argmax = make_shared<MultiVariatePoint<int>>(m_d,0.0);
    m_curentNeighbours.push_back(argmax);


    while (!m_curentNeighbours.empty() && iteration < m_maxIteration)
    {
        allAlphaLowerThanE = true;
        max = 0;
        if (debug)
        {
            Utils::separateur();
            displayPath();
            displayCurentNeighbours();
        }
        for (MultiVariatePointPtr<int> nu : m_curentNeighbours)
        {
            map<MultiVariatePointPtr<int>,double>::iterator it = m_alphaMap.find(nu);
            if (it != m_alphaMap.end())
                val = m_alphaMap[nu];
            else
                val =  computeLastAlphaNu(nu);
            if (val>e) allAlphaLowerThanE = false;
            nu->setAlpha(val);
            if (debug) cout << *nu << " " << val << " | ";
        }

        argmax = *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),(allAlphaLowerThanE) ? ageLess : alphaLess);
        max = argmax->getAlpha();

        m_path.push_back(argmax);
        m_alphaMap.insert(pair<MultiVariatePointPtr<int>,double>(argmax,max));
        addInterpolationPoint(getPoint(*argmax));

        if (m_method)
            updateNextPoints(argmax);
        else
            updateCurentNeighbours(argmax);

        iteration++;
        if (debug) cout << endl;

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI
        /*
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
        */
    }
    return iteration;
}
void Interpolation::updateCurentNeighbours(MultiVariatePointPtr<int> nu)
{
  m_curentNeighbours.remove(nu);
  for (MultiVariatePointPtr<int> mu : m_curentNeighbours)
      mu->incrWaitingTime();

  for (int k=0; k<nu->getD(); k++)
  {
      MultiVariatePointPtr<int> nu1 = make_shared<MultiVariatePoint<int>>(*nu);
      (*nu1)(k)++;
      if (isCorrectNeighbourToCurentPath(nu1))
          m_curentNeighbours.push_back(nu1);
  }
}
void Interpolation::updateNextPoints(MultiVariatePointPtr<int> nu)
{
    m_curentNeighbours.remove(nu);
    for (MultiVariatePointPtr<int> mu : m_curentNeighbours)
        mu->incrWaitingTime();

    vector<double> childrenVal;
    for (int i=0; i<m_d; i++)
    {
        childrenVal = Dichotomy::computeChildrenValue((*nu)(i));
        for (int k=0; k<int(childrenVal.size()); k++)
        {
            MultiVariatePointPtr<int> mu = make_shared<MultiVariatePoint<int>>(*nu);
            (*mu)(i) = Dichotomy::getIndice(childrenVal[k]);
            if (!indiceInPath(*mu) && !indiceInNeighborhood(*mu))
                m_curentNeighbours.push_back(mu);
        }
    }
}
bool Interpolation::isCorrectNeighbourToCurentPath(MultiVariatePointPtr<int> nu)
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
    for (MultiVariatePointPtr<int> nu : m_curentNeighbours)
        if (index==*nu)
            return true;
    return false;
}
bool Interpolation::indiceInPath(MultiVariatePoint<int> index)
{
    for (MultiVariatePointPtr<int> nu : m_path)
        if (index==*nu)
            return true;
    return false;
}
double Interpolation::computeLastAlphaNu(MultiVariatePointPtr<int> nu)
{
    double basisFuncProd;
    double res = Utils::gNd(getPoint(*nu));
    for (MultiVariatePointPtr<int> l : m_path)
    {
        if (l->lowerThan(*nu, m_method))
        {
            basisFuncProd = 1;
            for (int p=0; p<m_d; p++)
            {
                if (m_method == 0)
                    basisFuncProd *= lagrangeBasisFunction_1D((*l)(p),getPoint(*nu)(p),p);
                else if (m_method == 1)
                    basisFuncProd *= piecewiseFunction_1D((*l)(p),getPoint(*nu)(p),p); // TODO: getPoint(*) for method 1
            }
            res -= m_alphaMap[l] * basisFuncProd;
        }
    }
    m_alphaMap.insert(pair<MultiVariatePointPtr<int>,double>(nu,res));
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
        cout << " --> alpha" << *m_path[i] << " = " << m_alphaMap[m_path[i]] << endl;
    }
    cout << endl;
}
void Interpolation::displayAlphaTab()
{
    map<MultiVariatePointPtr<int>,double>::iterator it;
    cout << "Values of alpha (" << m_alphaMap.size() << ") :" << endl;
    for (it=m_alphaMap.begin(); it!=m_alphaMap.end(); it++)
        cout << "(" << *(get<0>(*it)) << ":" << getPoint(*(get<0>(*it))) << ":" << get<0>(*it) << ") : " << get<1>(*it) <<  endl;
}
void Interpolation::displayCurentNeighbours()
{
    cout << "Curent neighbours (" << m_curentNeighbours.size() << ") = ";
    for (MultiVariatePointPtr<int> nu : m_curentNeighbours)
        cout << "(" << (*nu) << ":" << getPoint(*nu) << ":" << nu << ") [" << nu->getWaitingTime() << "] | ";
    cout << endl << endl;
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
      for (MultiVariatePointPtr<int> nu : m_path)
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
          MultiVariatePoint<double> p;
          for (int j=0; j<nbPoints; j++)
          {
              file << x[j];
              for (MultiVariatePointPtr<int> nu : m_path)
              {
                  if (m_method == 0)
                      file << " " << lagrangeBasisFunction_1D((*nu)(0),x[j],0);
                  else if (m_method == 1)
                      file << " " << piecewiseFunction_1D((*nu)(0),x[j],0);
              }
              p = MultiVariatePoint<double>::toMonoVariatePoint(x[j]);
              file << " " << interpolation_ND(p);
              file << " " <<  Utils::gNd(p);
              file << endl;
          }
          for (MultiVariatePointPtr<int> nu : m_path)
              file << m_alphaMap[nu] << " ";
          file << endl;
          for (MultiVariatePoint<double> nu : m_interpolationNodes)
              file << nu(0) << " ";
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
        computeBoundaries(tk,&sup,&inf,axis);
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
double Interpolation::interpolation_ND(MultiVariatePoint<double>& x)
{
    double l_prod, sum = 0;
    for (MultiVariatePointPtr<int> nu : m_path)
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
