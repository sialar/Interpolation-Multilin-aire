#include "../include/Interpolation.hpp"

Interpolation::~Interpolation()
{}
Interpolation::Interpolation(vector<int> sizes)
{
    m_d = sizes.size();
    m_points.resize(m_d);
    m_middles.resize(m_d);
    for (int i=0; i<m_d; i++)
    {
        m_points[i].resize(sizes[i]);
        m_middles[i] = new BinaryTree();
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
    m_points[i] = Utils::createLejaSequence(nbPoints);
}

/******************************************************************************/


/************************* AI algo ********************************************/
int Interpolation::buildPathWithAIAlgo(int maxIteration, auto start_time, double threshold)
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
        updateCurentNeighbours(argmax);
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
            basisFuncProd *= lagrangeBasisFunction_1D(l(p),l(p),m_points[p][nu(p)],p);
        res -= m_alphaMap[l] * basisFuncProd;
    }
    m_alphaMap.insert(pair<MultiVariatePoint<int>,double>(nu,res));
    return res;
}
/******************************************************************************/


/*********************** Test functions ***************************************/
void Interpolation::testPathBuilt(int maxIteration, double threshold)
{
  auto start_time = chrono::steady_clock::now();
  int nbIterations = buildPathWithAIAlgo(maxIteration, start_time, threshold);
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
  for (int i=0; i<int(m_path.size()); i++)
  {
    cout << "(";
    for (int k=0; k<m_d-1; k++)
    cout << m_path[i](k) << ",";
    cout << m_path[i](m_d-1) << ")";
    if (i != int(m_path.size()-1))
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
        cout << nu << " ";
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
    double sup, inf, tk = m_points[axis][k];
    m_middles[axis]->searchNode(tk,&sup,&inf);
    if (t <= tk && t >= inf) return ((t-inf)/(tk-inf));
    else if (t >= tk && t <= sup) return ((t-sup)/(tk-sup));
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
            l_prod *= lagrangeBasisFunction_1D(nu(i),nu(i),x(i),i);
        sum += l_prod * m_alphaMap[nu];
    }
    return sum;
}
/******************************************************************************/
