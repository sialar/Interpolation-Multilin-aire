#include "../include/LagrangeInterpolation.hpp"

LagrangeInterpolation::LagrangeInterpolation(FunctionsPtr f, int nIter) :
    Interpolation(f,nIter)
{
    m_lejaSequence = Utils::loadLejaSequenceFromFile(m_maxIteration);
}

/************************* Data Points ****************************************/
MultiVariatePoint<double> LagrangeInterpolation::getPoint(MultiVariatePointPtr<int> nu)
{
	// Le point d'ordre nu est le point (m_lejaSequence[nu(0)], .., m_lejaSequence[nu(m_d)])
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

// Comparer des points multivariés selon les normes de leurs alpha (Erreur d'interpolation courante au point considéré)
bool alphaLess(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    return Utils::norm(nu->getAlpha(),2) < Utils::norm(mu->getAlpha(),2);
}
// Comparer des points multivariés selon le temps d'attente dans la liste des voisins courants (nombre d'itération passé dans la liste m_curentNeighbours)
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
 	// nu est le nouveau point d'interpolation (ce n'est plus un candidat)
	m_curentNeighbours.remove(nu);

	// Incrémenter le temps d'attente de tous les voisins (candidats)
	for (MultiVariatePointPtr<int> mu : m_curentNeighbours)
    	mu->incrWaitingTime();

	// Ajouter les nouveaus candidats (voisins de nu)
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
			// Si !indiceInPath(mu), nu ne peux pas être un candidat (il faut concerver un ensemble monotone)
            if (!indiceInPath(mu)) return false;
            mu(i)++;
        }
    }
    return true;
}

// Vérifier si index est déjà un point d'interpolation (ajouté dans path)
bool LagrangeInterpolation::indiceInPath(MultiVariatePoint<int> index)
{
    for (MultiVariatePointPtr<int> nu : m_path)
        if (index==*nu)
            return true;
    return false;
}
/******************************************************************************/

/*********************** Interpolation ****************************************/
// Polynôme hiérarchique de Lagrange
double LagrangeInterpolation::basisFunction_1D(int k, double t, int axis)
{
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        prod *= (t-m_lejaSequence[i]) / (m_lejaSequence[k]-m_lejaSequence[i]) ;
    return prod;
}
/******************************************************************************/
