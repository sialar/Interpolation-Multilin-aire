#include "../include/LagrangeInterpolation.hpp"

LagrangeInterpolation::LagrangeInterpolation(FunctionsPtr f, int nIter) :
    Interpolation(f,nIter)
{
    // Charger la séquence de points de Leja 1d depuis le fichier data/leja_sequence.dat
    m_lejaSequence = Utils::loadLejaSequenceFromFile(m_maxIteration);
}


/******************************************************************************/
/************************ Points d'interpolation ******************************/
MultiVariatePoint<double> LagrangeInterpolation::getPoint(MultiVariatePointPtr<int> nu)
{
	 // Le point d'ordre nu est le point (m_lejaSequence[nu(0)], .., m_lejaSequence[nu(m_d)])
    MultiVariatePoint<double> point(m_d,0,0.0);
    for (int i=0; i<m_d; i++)
        // sur chaque direction i, on prend le point ième point 1d de la séquence de Leja
        point(i) = m_lejaSequence[(*nu)(i)];
    return point;
}
void LagrangeInterpolation::addInterpolationPoint(MultiVariatePoint<double>p)
{
    // Ajout dans m_interpolationNodes
    m_interpolationNodes.push_back(p);

    // Ajout de chaque coordonnée i dans m_interpolationPoints[i]
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
/*************************** Algorithme AI ************************************/
bool alphaLess(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    // Comparer des points multivariés selon les normes de leurs alpha (Erreur d'interpolation courante au point considéré)
    return Utils::norm(nu->getAlpha(),2) < Utils::norm(mu->getAlpha(),2);
}
bool ageLess(MultiVariatePointPtr<int> nu, MultiVariatePointPtr<int> mu)
{
    // Comparer des points multivariés selon le temps d'attente dans la liste des voisins courants (nombre d'itération passé dans la liste m_curentNeighbours)
    return nu->getNbWaitingIter() < mu->getNbWaitingIter();
}
MultiVariatePointPtr<int> LagrangeInterpolation::maxElement(int iteration)
{
  // 1 fois sur 2, on choisi le point où l'erreur est la plus élevé
  if (iteration%2)
      return *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),alphaLess);
  // 1 fois sur 2, on choisi le point qui a le plus attendu dans la liste des voisins
  // si il y en a plusieurs, on en choisit un où l'erreur est nulle pour éviter de bloquer une direction
  // En effet, si l'erreur en un point p vaut 0 elle le restaura si on ne procede que par adaptivité
  // Or parfois, une erreur peut être nulle par hasard (dépend de l'interpolé)
  // Exemple f(x,y) = sqrt(1-x**2)
  else
  {
      MultiVariatePointPtr<int> mu = *max_element(m_curentNeighbours.begin(), \
                                          m_curentNeighbours.end(),ageLess);
      for (MultiVariatePointPtr<int> nu : m_curentNeighbours)
          if (nu->alphaIsNull() && nu->getNbWaitingIter()==mu->getNbWaitingIter())
              return nu;
      return mu;
  }
}
MultiVariatePointPtr<int> LagrangeInterpolation::getFirstMultivariatePoint()
{
    // Le plus petit indice est 0
    MultiVariatePointPtr<int> nu = make_shared<MultiVariatePoint<int>>(m_d,m_n,0);
    return nu;
}
void LagrangeInterpolation::updateCurentNeighbours(MultiVariatePointPtr<int> nu)
{
 	// nu est le nouveau point d'interpolation (ce n'est plus un candidat)
	m_curentNeighbours.remove(nu);

	// Incrémenter le temps d'attente de tous les voisins (candidats)
	for (MultiVariatePointPtr<int> mu : m_curentNeighbours)
    	mu->incrNbWaitingIter();

	// Ajouter les nouveaus candidats (voisins de nu)
  	for (int k=0; k<m_d; k++)
  	{
    	MultiVariatePointPtr<int> nu1 = make_shared<MultiVariatePoint<int>>(*nu);
      // nu1 est un voisin de nu si pour tout i nu1(i) = nu(i) sauf pour un i0 où nu1(i0) = nu(i0)+1
    	(*nu1)(k)++;
      // Vérifier que nu1 est un voisin de m_path courant
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
bool LagrangeInterpolation::indiceInPath(MultiVariatePoint<int> index)
{
    // Vérifier si index est déjà un point d'interpolation (ajouté dans path)
    for (MultiVariatePointPtr<int> nu : m_path)
        if (index==*nu)
            return true;
    return false;
}

/******************************************************************************/
/************************** Fonctions de base *********************************/
double LagrangeInterpolation::basisFunction_1D(int k, double t, int axis)
{
    // Polynôme hiérarchique de Lagrange (voir formule h_k page 10 du rapport)
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        prod *= (t-m_lejaSequence[i]) / (m_lejaSequence[k]-m_lejaSequence[i]) ;
    return prod;
}
