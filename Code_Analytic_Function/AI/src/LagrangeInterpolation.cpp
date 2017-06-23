#include "../include/LagrangeInterpolation.hpp"

LagrangeInterpolation::LagrangeInterpolation(int d, int n, int nIter, Function f) :
    Interpolation(d,n,nIter,f)
{
    m_lejaSequence = Utils::loadLejaSequenceFromFile(m_maxIteration);
}

/************************* Data Points ****************************************/
MultiVariatePoint<double> LagrangeInterpolation::getPoint(MultiVariatePointPtr<int> nu)
{
    MultiVariatePoint<double> point(m_d,0,0.0);
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
bool alphaLess(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    return Utils::norm(nu->getAlpha(),2) < Utils::norm(mu->getAlpha(),2);
}
bool ageLess(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    return nu->getWaitingTime() < mu->getWaitingTime();
}
MultiVariatePointPtr<int> LagrangeInterpolation::maxElement(int iteration)
{
  if (iteration%4)
      return *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),alphaLess);
  else
  {
      MultiVariatePointPtr<int> mu = *max_element(m_curentNeighbours.begin(), \
                                          m_curentNeighbours.end(),ageLess);
      for (MultiVariatePointPtr<int> nu : m_curentNeighbours)
          if (nu->alphaIsNull() && nu->getWaitingTime()==mu->getWaitingTime())
              return nu;
      return mu;
  }
}
MultiVariatePointPtr<int> LagrangeInterpolation::getFirstMultivariatePoint()
{
    MultiVariatePointPtr<int> nu = make_shared<MultiVariatePoint<int>>(m_d,m_n,0);
    return nu;
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
    MultiVariatePoint<int> mu(*nu);
    string temp;
    for (int i=0; i<m_d; i++)
    {
        if ((*nu)(i))
        {
            mu(i)--;
            if (!indiceInPath(mu)) return false;
            mu(i)++;
        }
    }
    return true;
}
bool LagrangeInterpolation::indiceInPath(MultiVariatePoint<int> index)
{
    for (MultiVariatePointPtr<int> nu : m_path)
        if (index==*nu)
            return true;
    return false;
}
/******************************************************************************/

/************************* Display functions **********************************/
void LagrangeInterpolation::saveInterpolationBasisFunctions()
{
  ofstream file(Utils::projectPath + "data/basis_functions.txt", ios::out | ios::trunc);
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
            file << " " << basisFunction_1D((*nu)(0),x[j],0);
        p = MultiVariatePoint<double>::toMonoVariatePoint(x[j]);
        file << " " <<  func(p)[0];
        file << endl;
      }
      for (MultiVariatePointPtr<int> nu : m_path)
          file << nu->getAlpha()[0] << " ";
      file << endl;
      for (MultiVariatePoint<double> nu : m_interpolationNodes)
      file << nu(0) << " ";
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}
/******************************************************************************/


/*********************** Interpolation ****************************************/
double LagrangeInterpolation::basisFunction_1D(int k, double t, int axis)
{
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        prod *= (t-m_lejaSequence[i]) / (m_lejaSequence[k]-m_lejaSequence[i]) ;
    return prod;
}
/******************************************************************************/