#include "../include/LagrangeInterpolation.hpp"

LagrangeInterpolation::LagrangeInterpolation(int d, int nIter) :
    Interpolation(d,nIter)
{
    m_lejaSequence = Utils::loadLejaSequenceFromFile(m_maxIteration);
}

/************************* Data Points ****************************************/
MultiVariatePoint<double> LagrangeInterpolation::getPoint(MultiVariatePointPtr<int> nu)
{
    MultiVariatePoint<double> point(m_d,0.0);
    for (int i=0; i<m_d; i++)
        point(i) = m_lejaSequence[(*nu)(i)];
    return point;
}
void LagrangeInterpolation::addInterpolationPoint(MultiVariatePoint<double>p)
{
    m_interpolationNodes.push_back(p);

    bool found;
    for (int i=0; i<m_d; i++)
    {
        found = false;
        for (double x : m_interpolationPoints[i])
            if (x == p(i)) found = true;
        if (!found)
            m_interpolationPoints[i].push_back(p(i));
    }
}
/******************************************************************************/


/************************* AI algo ********************************************/
bool alphaLess1(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    return abs(nu->getAlpha()) < abs(mu->getAlpha());
}
bool ageLess1(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    return nu->getWaitingTime() < mu->getWaitingTime();
}
int LagrangeInterpolation::buildPathWithAIAlgo(auto start_time, double threshold, bool debug)
{
    m_curentNeighbours.clear();
    m_path.clear();
    MultiVariatePointPtr<int> argmax = make_shared<MultiVariatePoint<int>>(m_d,0.0);
    m_curentNeighbours.push_back(argmax);
    int iteration = 0;
    double val;

    while (!m_curentNeighbours.empty() && iteration < m_maxIteration)
    {
        if (debug)
        {
            Utils::separateur();
            displayPath();
            displayCurentNeighbours();
        }
        for (MultiVariatePointPtr<int> nu : m_curentNeighbours)
        {
            if (nu->alphaAlreadyComputed()) val = nu->getAlpha();
            else val =  computeLastAlphaNu(nu);
            if (debug) cout << *nu << " " << val << " | ";
        }
        if (debug) cout << endl;
        argmax = *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),(iteration%4) ? alphaLess1 : ageLess1);
        m_path.push_back(argmax);
        addInterpolationPoint(getPoint(argmax));
        updateCurentNeighbours(argmax);
        iteration++;

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI
        if (m_saveError)
        {
            auto end_time = chrono::steady_clock::now();
            std::chrono::duration<double> run_time = end_time - start_time;
            double error = tryWithCurentPath();
            m_errors.insert(pair<int, double>(iteration, tryWithCurentPath()));
            if (iteration%(m_maxIteration/10)==0)
            {
                cout << "   - Interpolation error after " << iteration << " iterations: " << error;
                cout << " | Elapsed time : "  << run_time.count() << endl;
            }
            if (error < threshold)
            {
                cout << endl << "   - AI Algo stop after " << iteration << " iterations";
                cout << " | Elapsed time : "  << run_time.count() << endl;
                return iteration;
            }
            saveErrorsInFile();
        }
    }
    return iteration;
}
void LagrangeInterpolation::updateCurentNeighbours(MultiVariatePointPtr<int> nu)
{
  m_curentNeighbours.remove(nu);
  for (MultiVariatePointPtr<int> mu : m_curentNeighbours)
      mu->incrWaitingTime();

  for (int k=0; k<m_d; k++)
  {
      MultiVariatePointPtr<int> nu1 = make_shared<MultiVariatePoint<int>>(*nu);
      (*nu1)(k)++;
      if (isCorrectNeighbourToCurentPath(nu1))
          m_curentNeighbours.push_back(nu1);
  }
}
bool LagrangeInterpolation::isCorrectNeighbourToCurentPath(MultiVariatePointPtr<int> nu)
{
    bool isNeighbour[m_d];
    for (int i=0; i<m_d; i++)
        isNeighbour[i] = true;

    MultiVariatePoint<int> index(*nu);
    for (int i=0; i<m_d; i++)
        if ((*nu)(i))
        {
            index(i)--;
            isNeighbour[i] = indiceInPath(index) && !indiceInNeighborhood(index);
            index(i)++;
        }
    bool res = true;
    for (int i=0; i<m_d; i++)
        res = res && isNeighbour[i];
    return res;
}
bool LagrangeInterpolation::indiceInNeighborhood(MultiVariatePoint<int> index)
{
    for (MultiVariatePointPtr<int> nu : m_curentNeighbours)
        if (index==*nu)
            return true;
    return false;
}
bool LagrangeInterpolation::indiceInPath(MultiVariatePoint<int> index)
{
    for (MultiVariatePointPtr<int> nu : m_path)
        if (index==*nu)
            return true;
    return false;
}
double LagrangeInterpolation::computeLastAlphaNu(MultiVariatePointPtr<int> nu)
{
    double basisFuncProd;
    double res = Utils::gNd(getPoint(nu));
    for (MultiVariatePointPtr<int> l : m_path)
    {
        basisFuncProd = 1;
        for (int p=0; p<m_d; p++)
            basisFuncProd *= lagrangeBasisFunction_1D((*l)(p),getPoint(nu)(p),p);
        res -= l->getAlpha() * basisFuncProd;
    }
    m_alphaMap.insert(pair<MultiVariatePointPtr<int>,double>(nu,res));
    nu->setAlpha(res);
    return res;
}
/******************************************************************************/


/*********************** Test functions ***************************************/
void LagrangeInterpolation::testPathBuilt(double threshold, bool debug)
{
  auto start_time = chrono::steady_clock::now();
  int nbIterations = buildPathWithAIAlgo(start_time, threshold, debug);
  auto end_time = chrono::steady_clock::now();
  std::chrono::duration<double> run_time = end_time - start_time;
  cout << endl << " - Time required to compute the path with " << nbIterations <<
  " iterations: " << run_time.count() << "s" << endl;
}
double LagrangeInterpolation::tryWithCurentPath()
{
  vector<double> realValues, estimate;
  for (MultiVariatePoint<double> p : m_testPoints)
  {
    realValues.push_back(Utils::gNd(p));
    estimate.push_back(interpolation_ND(p,m_path.size()));
  }
  return Utils::interpolationError(realValues,estimate);
}
/******************************************************************************/


/************************* Display functions **********************************/
void LagrangeInterpolation::displayPath()
{
    // Format: (nu:points[nu]) --> () ...
    cout << "Chemin =";
    int n = m_path.size();
    for (int i=0; i<n; i++)
    {
        if (i>0) cout << "\t";
        cout << " " << i << " :";
        cout << " [" << *m_path[i] << ":" << getPoint(m_path[i]) << ":" << m_path[i] << "]";
        cout << " --> alpha" << *m_path[i] << " = " << m_path[i]->getAlpha() << endl;
    }
    cout << endl;
}
void LagrangeInterpolation::displayAlphaTab()
{
    map<MultiVariatePointPtr<int>,double>::iterator it;
    cout << "Values of alpha (" << m_alphaMap.size() << ") :" << endl;
    for (it=m_alphaMap.begin(); it!=m_alphaMap.end(); it++)
        cout << "(" << *(get<0>(*it)) << ":" << getPoint((get<0>(*it))) << ":" << get<0>(*it) << ") : " << get<1>(*it) <<  endl;
}
void LagrangeInterpolation::displayCurentNeighbours()
{
    cout << "Curent neighbours (" << m_curentNeighbours.size() << ") = ";
    for (MultiVariatePointPtr<int> nu : m_curentNeighbours)
        cout << "(" << (*nu) << ":" << getPoint(nu) << ":" << nu << ") [" << nu->getWaitingTime() << "] | ";
    cout << endl << endl;
}

void LagrangeInterpolation::savePathInFile()
{
  ofstream file("data/path.txt", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==2)
    {
      //cout << " - The path is saved in data/path.txt" << endl;
      file << m_interpolationPoints[0].size() << " " <<  m_interpolationPoints[1].size() << endl;
      file << m_path.size() << endl;
      MultiVariatePoint<double> x(m_d,0.0);
      for (MultiVariatePointPtr<int> nu : m_path)
      {
        x = getPoint(nu);
        file << x(0) << " " << x(1) << " ";
        file << (*nu)(0) << " " << (*nu)(1) << endl;
      }
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}
void LagrangeInterpolation::storeInterpolationBasisFunctions()
{
  ofstream file("data/basis_functions.txt", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==1)
    {
      vector<double> x;
      int nbPoints = int(m_testPoints.size());
      for (int i=0; i<nbPoints; i++)
          x.push_back(m_testPoints[i](0));
      sort(x.begin(), x.end());
      file << m_path.size() << " " << nbPoints << " " << "0" << endl;
      MultiVariatePoint<double> p;
      for (int j=0; j<nbPoints; j++)
      {
        file << x[j];
        for (MultiVariatePointPtr<int> nu : m_path)
            file << " " << lagrangeBasisFunction_1D((*nu)(0),x[j],0);
        p = MultiVariatePoint<double>::toMonoVariatePoint(x[j]);
        file << " " <<  Utils::gNd(p);
        file << endl;
      }
      for (MultiVariatePointPtr<int> nu : m_path)
      file << nu->getAlpha() << " ";
      file << endl;
      for (MultiVariatePoint<double> nu : m_interpolationNodes)
      file << nu(0) << " ";
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}

void LagrangeInterpolation::storeInterpolationProgression()
{
  ofstream file("data/interpolation_progression.txt", ios::out | ios::trunc);
  if(file)
  {
      if (m_d==1)
      {
          vector<double> x;
          int nbPoints = int(m_testPoints.size());
          for (int i=0; i<nbPoints; i++)
              x.push_back(m_testPoints[i](0));
          sort(x.begin(), x.end());
          MultiVariatePoint<double> p;
          vector<double> tempPath;
          for (double t : x)
          {
              for (int i=0; i<int(m_path.size()); i++)
              {
                  p = MultiVariatePoint<double>::toMonoVariatePoint(t);
                  file << interpolation_ND(p, i+1) << " ";
              }
              file << endl;
          }
      }
      file.close();
  }
  else
      cerr << "Error while opening the file!" << endl;
}

/******************************************************************************/


/*********************** Interpolation ****************************************/
double LagrangeInterpolation::lagrangeBasisFunction_1D(int k, double t, int axis)
{
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        prod *= (t-m_lejaSequence[i]) / (m_lejaSequence[k]-m_lejaSequence[i]) ;
    return prod;
}
double LagrangeInterpolation::interpolation_ND(MultiVariatePoint<double>& x, int end)
{
  double l_prod, sum = 0;
  for (int k=0; k<end; k++)
  {
      l_prod = 1;
      for (int i=0; i<m_d; i++)
          l_prod *= lagrangeBasisFunction_1D((*m_path[k])(i),x(i),i);
      sum += l_prod * m_path[k]->getAlpha();
  }
  return sum;
}
/******************************************************************************/
