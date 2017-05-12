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
        {
            m_middles[i] = new BinaryTree();
            m_middles[i]->initTree(int(log(sizes[i])/log(2))-1);
        }
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
      m_points[i].push_back(-1);
      m_points[i].push_back(1);
      m_points[i].push_back(0);
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

/******************************************************************************/


/************************* AI algo ********************************************/
int Interpolation::buildPathWithAIAlgo(int maxIteration, auto start_time, double threshold, bool debug)
{
    m_curentNeighbours.clear();
    m_path.clear();

    MultiVariatePoint<int> nu(m_d,0);
    m_path.push_back(nu);
    m_alphaMap[nu] = Utils::gNd(getPoint(nu));

    double val, max;
    MultiVariatePoint<int> argmax(m_d,0);
    int iteration = 1;

    while (iteration < maxIteration)
    {
        max = -numeric_limits<double>::max();
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


            if (val >= max)
            {
                max = val;
                argmax = nu;
            }
        }
        m_path.push_back(argmax);
        iteration++;

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI
        if (iteration%(maxIteration/10)==0)
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
        if (nu(i) == 0)
        {
            mu(i) = 1;
            if (!indiceInPath(mu))
                m_curentNeighbours.push_back(mu);
        }
        else if (nu(i) == 1)
        {
            MultiVariatePoint<int> mu(nu);
            mu(i) = 2;
            if (!indiceInPath(mu))
                m_curentNeighbours.push_back(mu);
        }
        else
        {
            double inf, sup;
            Node* node = m_middles[i]->searchNode(getPoint(nu)(i),&inf,&sup);
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
    double res = Utils::gNd(getPoint(nu));
    for (MultiVariatePoint<int> l : m_path)
    {
        basisFuncProd = 1;
        for (int p=0; p<m_d; p++)
        {
            if (m_method == 0)
                basisFuncProd *= lagrangeBasisFunction_1D(l(p),l(p),m_points[p][nu(p)],p);
            else if (m_method == 1)
                basisFuncProd *= piecewiseFunction_1D(l(p),m_points[p][nu(p)],p);
            else
                basisFuncProd *= piecewiseLagrangeBasisFunction_1D(l(p),l(p),m_points[p][nu(p)],p);
        }
        res -= m_alphaMap[l] * basisFuncProd;
    }
    m_alphaMap.insert(pair<MultiVariatePoint<int>,double>(nu,res));
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
  cout << "Chemin = ";
  int n = min(10,int(m_path.size()));
  for (int i=0; i<n; i++)
  {
    cout << "(";
    for (int k=0; k<m_d-1; k++)
        cout << " [" << m_path[i](k) << ",";
    cout << m_path[i](m_d-1) << "] = ";
    for (int k=0; k<m_d-1; k++)
        cout << "(" << getPoint(m_path[i])(k) << ",";
    cout << getPoint(m_path[i])(m_d-1) << ") )";
    if (i != n-1)
    cout << " --> ";
  }
  cout << endl;
}
void Interpolation::displayAlphaTab()
{
    map<MultiVariatePoint<int>,double>::iterator it;
    cout << "Values of alpha (" << m_alphaMap.size() << ") :" << endl;
    for (it=m_alphaMap.begin(); it!=m_alphaMap.end(); it++)
        cout << get<0>(*it) << ": " << get<1>(*it) << endl;
}
void Interpolation::displayCurentNeighbours()
{
    cout << "Curent neighbours (" << m_curentNeighbours.size() << ") = ";
    for (MultiVariatePoint<int> nu : m_curentNeighbours)
        cout << getPoint(nu) << " ";
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
double Interpolation::piecewiseFunction_1D(int k, double t, int axis)
{
    if (k==0) return 1;
    if (k==1) return -0.5*(t+1);
    double sup, inf, tk = m_points[axis][k];
    m_middles[axis]->searchNode(tk,&sup,&inf);
    if (t <= tk && t >= inf) return ((t-inf)/(tk-inf));
    else if (t >= tk && t <= sup) return ((t-sup)/(tk-sup));
    else return 0;
}
double Interpolation::piecewiseLagrangeBasisFunction_1D(int j, int k, double t, int axis)
{
    if (k==0) return 1;
    if (k==1) return -0.25*(t+1)*(t-3);
    double sup, inf, tj = m_points[axis][j];
    m_middles[axis]->searchNode(tj,&sup,&inf);
    if (t <= sup && t >= inf)
        return -4 * (t-inf) * (t-sup) / pow(sup-inf,2);
    else return 0;
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
            else
                l_prod *= piecewiseLagrangeBasisFunction_1D(nu(i),nu(i),x(i),i);
        }
        sum += l_prod * m_alphaMap[nu];
    }
    return sum;
}
/******************************************************************************/
