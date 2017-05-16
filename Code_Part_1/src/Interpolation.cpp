#include "../include/Interpolation.hpp"

Interpolation::~Interpolation()
{}
Interpolation::Interpolation(vector<int> sizes, int method)
{
    m_method = method;
    m_d = sizes.size();
    m_points.resize(m_d);
    m_middles.resize(m_d);
    for (int i=0; i<m_d; i++)
    {
        if (m_method)
            m_middles[i] = new BinaryTree(int(log(sizes[i])/log(2)));
        else
            m_points[i].resize(sizes[i]);
    }
}

/************************* Data Points ****************************************/
MultiVariatePoint<double> Interpolation::getPoint(MultiVariatePoint<int> nu)
{
    MultiVariatePoint<double> point(m_d,0.0);
    for (int i=0; i<m_d; i++)
        point(i) = m_points[i][nu(i)];
    return point;
}
void Interpolation::setDirPoints(int i, int nbPoints)
{
    if (m_method)
    {
      if (nbPoints > 0) m_points[i].push_back(0);
      if (nbPoints > 1) m_points[i].push_back(-1);
      if (nbPoints > 2) m_points[i].push_back(1);
      int e = 0, k = 3;
      double denom;
      while (k<nbPoints)
      {
          e++;
          denom = pow(2,e);
          for (int j=denom-1; j>=0 && int(m_points[i].size())<nbPoints; j--)
              if (j%2) m_points[i].push_back(-j/denom);
          for (int j=0; j<denom && int(m_points[i].size())<nbPoints; j++)
              if (j%2) m_points[i].push_back(j/denom);
          k += denom;
      }
    }
    else
        m_points[i] = Utils::createLejaSequence(nbPoints);
}
MultiVariatePoint<double> Interpolation::getMultiPoint(MultiVariatePoint<int> nu)
{
    MultiVariatePoint<double> x(nu.getD(),0.0);
    for (int i=0; i<nu.getD(); i++)
        x(i) = BinaryTree::getValue(nu(i));
    return x;
}
MultiVariatePoint<int> Interpolation::getMultiIndice(MultiVariatePoint<double> x)
{
    MultiVariatePoint<int> nu(x.getD(),0);
    for (int i=0; i<x.getD(); i++)
        nu(i) = BinaryTree::getIndice(x(i));
    return nu;
}
/******************************************************************************/


/************************* AI algo ********************************************/
int Interpolation::buildPathWithAIAlgo(int maxIteration, auto start_time, double threshold, bool debug)
{
    m_curentNeighbours.clear();
    m_path.clear();

    MultiVariatePoint<int> nu(m_d,0);
    m_path.push_back(nu);
    m_alphaMap.insert(pair<MultiVariatePoint<int>,double>(nu,Utils::gNd(getPoint(nu))));

    double val, max;
    MultiVariatePoint<int> argmax(m_d,0);
    m_curentNeighbours.push_back(argmax);
    int iteration = 1;

    while (!m_curentNeighbours.empty() && iteration < maxIteration)
    {
        max = 0;
        if (m_method)
            updateNextPoints(argmax);
        else
            updateCurentNeighbours(argmax);

        if (debug)
        {
            displayPath();
            displayCurentNeighbours();
        }

        for (MultiVariatePoint<int> nu : m_curentNeighbours)
        {
            map<MultiVariatePoint<int>,double>::iterator it = m_alphaMap.find(nu);
            if (it != m_alphaMap.end())
                val = m_alphaMap[nu];
            else
                val =  computeLastAlphaNu(nu);

            if (debug) cout << nu << " " << val << " | ";
            if (abs(val) >= abs(max))
            {
                max = val;
                argmax = nu;
            }
        }
        if (debug) cout << endl;

        m_path.push_back(argmax);
        m_alphaMap.insert(pair<MultiVariatePoint<int>,double>(argmax,max));
        iteration++;

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI

        if ((maxIteration>10) && iteration%(maxIteration/10)==0)
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
    return iteration;
}
void Interpolation::updateCurentNeighbours(MultiVariatePoint<int>& nu)
{
  m_curentNeighbours.remove(nu);
  MultiVariatePoint<int> nu1(nu);
  int max_indice[m_d];
  for (int i=0; i<m_d; i++)
  max_indice[i] = m_points[i].size()-1;
  for (int k=0; k<nu.getD(); k++)
  {
    if (nu1(k)<max_indice[k])
    {
      nu1(k)++;
      if (isCorrectNeighbourToCurentPath(nu1))
      m_curentNeighbours.push_back(nu1);
      nu1(k)--;
    }
  }
}
void Interpolation::updateNextPoints(MultiVariatePoint<int>& nu)
{
    m_curentNeighbours.remove(nu);

    for (int i=0; i<m_d; i++)
    {
        MultiVariatePoint<int> mu(nu);
        double inf, sup;
        Node* node = m_middles[i]->searchNode(getPoint(nu)(i),&inf,&sup,true);
        int inf_index, sup_index;
        if (node->left())
        {
            inf = node->left()->key();
            inf_index = BinaryTree::getIndice(inf);
            mu(i) = inf_index;
            if (!indiceInPath(mu))
                m_curentNeighbours.push_back(mu);
        }
        if (node->right())
        {
            sup = node->right()->key();
            sup_index = BinaryTree::getIndice(sup);
            mu(i) = sup_index;
            if (!indiceInPath(mu))
                m_curentNeighbours.push_back(mu);
        }
    }
}
bool Interpolation::isCorrectNeighbourToCurentPath(MultiVariatePoint<int>& nu)
{
    bool isNeighbour[m_d];
    for (int i=0; i<m_d; i++)
        isNeighbour[i] = true;

    MultiVariatePoint<int> index(nu);
    for (int i=0; i<m_d; i++)
        if (nu(i))
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
bool Interpolation::indiceInPath(MultiVariatePoint<int>& index)
{
    for (MultiVariatePoint<int> nu : m_path)
        if (index==nu)
            return true;
    return false;
}
double Interpolation::computeLastAlphaNu(MultiVariatePoint<int>& nu)
{
    double basisFuncProd;
    //cout << "\n\tCalcul de alpha" << nu << "" << endl;
    double res = Utils::gNd(getPoint(nu));
    //cout << "\tf" << nu << " = " << res << endl;
    //cout << "\tLongeur du chemin = " << m_path.size() << endl;

    for (MultiVariatePoint<int> l : m_path)
    {
        basisFuncProd = 1;
        for (int p=0; p<m_d; p++)
        {
            if (m_method == 0)
                basisFuncProd *= lagrangeBasisFunction_1D(l(p),l(p),m_points[p][nu(p)],p);
            else if (m_method == 1)
                basisFuncProd *= piecewiseFunction_1D(l(p),m_points[p][nu(p)],p);
        }
        //cout << "\t" << l << " " << m_alphaMap[l] << " " << basisFuncProd << endl;
        res -= m_alphaMap[l] * basisFuncProd;
        //cout << "\talpha"  << nu << "=" << res << endl;
    }
    return res;
}
/******************************************************************************/


/*********************** Test functions ***************************************/
void Interpolation::testPathBuilt(int maxIteration, double threshold, bool debug)
{
  auto start_time = chrono::steady_clock::now();
  int nbIterations = buildPathWithAIAlgo(maxIteration, start_time, threshold, debug);
  auto end_time = chrono::steady_clock::now();
  std::chrono::duration<double> run_time = end_time - start_time;
  cout << endl << " - Time required to compute the path with " << nbIterations <<
  " iterations: " << run_time.count() << "s" << endl;
}
double Interpolation::testAlphaNuComputation(MultiVariatePoint<int>& nu)
{
    auto start_time = chrono::steady_clock::now();
    double val = computeLastAlphaNu(nu);
    chrono::duration<double> run_time = chrono::steady_clock::now() - start_time;
    cout << " - Time required to compute alpha(" << nu << "): " << run_time.count() << " s";
    cout << " Correct Value = " << val << endl;
    return val;
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
        cout << " [" << m_path[i] << ":" << getPoint(m_path[i]) << "]";
        cout << " --> alpha(" << m_path[i] << ") = " << m_alphaMap[m_path[i]] << endl;
    }
}
void Interpolation::displayAlphaTab()
{
    map<MultiVariatePoint<int>,double>::iterator it;
    cout << "Values of alpha (" << m_alphaMap.size() << ") :" << endl;
    for (it=m_alphaMap.begin(); it!=m_alphaMap.end(); it++)
        cout << "(" << get<0>(*it) << ":" << getPoint(get<0>(*it)) << ") : " << get<1>(*it) <<  endl;
}
void Interpolation::displayCurentNeighbours()
{
    cout << "Curent neighbours (" << m_curentNeighbours.size() << ") = ";
    for (MultiVariatePoint<int> nu : m_curentNeighbours)
        cout << "(" << nu << ":" << getPoint(nu) << ") | ";
    cout << endl;
}
void Interpolation::displayInterpolationPoints()
{
    cout << " - Dimension d = " << m_d << endl;
    for (int i=0; i<m_d; ++i)
    {
        cout << " - " << m_points[i].size() << " points in direction " << i << " : { ";
        for (size_t j=0; j<m_points[i].size()-1; ++j)
            cout << m_points[i][j] << " ; ";
        cout << m_points[i][m_points[i].size()-1] << " }" << endl;
    }
}
void Interpolation::savePathInFile()
{
  ofstream file("python/output_algo_AI.txt", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==2)
    {
      cout << " - The path is saved in python/output_algo_AI.txt" << endl;
      file << m_points[0].size() << " " <<  m_points[1].size() << endl;
      file << m_path.size() << endl;
      for (MultiVariatePoint<int> nu : m_path)
      file << nu(0) << " " << nu(1) << endl;
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
              for (MultiVariatePoint<int> nu : m_path)
              {
                  if (m_method == 0)
                      file << " " << lagrangeBasisFunction_1D(nu(0),m_path.size(),x[j],0);
                  else if (m_method == 1)
                      file << " " << piecewiseFunction_1D(nu(0),x[j],0);
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
        double sup, inf, tk = m_points[axis][k];
        m_middles[axis]->searchNode(tk,&sup,&inf, false);
        if (t <= tk && t >= inf) return ((t-inf)/(tk-inf));
        else if (t >= tk && t <= sup) return ((t-sup)/(tk-sup));
        else return 0;
    }
}
double Interpolation::lagrangeBasisFunction_1D(int j, int k, double t, int axis)
{
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        if (i!=j)
            prod *= (t-m_points[axis][i]) / (m_points[axis][j]-m_points[axis][i]) ;
    return prod;
}
double Interpolation::interpolation_ND(MultiVariatePoint<double> x)
{
    double l_prod, sum = 0;
    for (MultiVariatePoint<int> nu : m_path)
    {
        l_prod = 1;
        for (int i=0; i<m_d; i++)
        {
            if (m_method == 0)
                l_prod *= lagrangeBasisFunction_1D(nu(i),nu(i),x(i),i);
            else if (m_method == 1)
                l_prod *= piecewiseFunction_1D(nu(i),x(i),i);
        }
        sum += l_prod * m_alphaMap[nu];
    }
    return sum;
}
/******************************************************************************/
