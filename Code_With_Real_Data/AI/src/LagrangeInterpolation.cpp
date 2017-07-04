#include "../include/LagrangeInterpolation.hpp"

LagrangeInterpolation::LagrangeInterpolation(int d, string core, vector<string> cs, int nIter) :
    Interpolation(d, core, cs ,nIter)
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
  if (iteration%2)
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
