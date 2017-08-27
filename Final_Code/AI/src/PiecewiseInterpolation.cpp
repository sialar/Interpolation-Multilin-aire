#include "../include/PiecewiseInterpolation.hpp"

PiecewiseInterpolation::PiecewiseInterpolation(int d, int n, int nIter, int method) :
    Interpolation(d, n, nIter), m_method(method)
{
    m_trees.resize(m_d);
    for (int i=0; i<m_d; i++)
      m_trees[i] = make_shared<BinaryTree>();
  }

void PiecewiseInterpolation::clearAllTrees()
{
    for (int i=0; i<int(m_trees.size()); i++)
        m_trees[i]->clearTree();
}

/************************* Data Points ****************************************/
// Le point d'ordre nu est le point (x[nu(0)], .., x[nu(m_d)]), nu[i] est un code de Huffman
MultiVariatePoint<double> PiecewiseInterpolation::getPoint(MultiVariatePointPtr<string> nu)
{
    MultiVariatePoint<double> point(m_d,0,0.0);
    for (int i=0; i<m_d; i++)
        point(i) = BinaryTree::getValueFromCode((*nu)(i));
    return point;
}
void PiecewiseInterpolation::addInterpolationPoint(MultiVariatePoint<double>p)
{
    m_interpolationNodes.push_back(p);

    bool found;
    for (int i=0; i<m_d; i++)
    {
        found = false;
        for (double x : m_interpolationPoints[i])
            if (x == p(i)) found = true;
        if (!found)
        {
            m_interpolationPoints[i].push_back(p(i));
            m_trees[i]->addNode(p(i));
        }
    }
}

// Recherche des valeurs (inf et sup) les plus proches de t (de part et d'autre de t) sur la direction axis
void PiecewiseInterpolation::computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis)
{
    /*
    vector<double> higherPoints, lowerPoints;
    for (double x : m_interpolationPoints[axis])
    {
        if (x>t) higherPoints.push_back(x);
        else if (x<t) lowerPoints.push_back(x);
        else break;
        // break when t is found, do not look at children when looking for boundaries
        // because children were constructed after t
    }
    if (higherPoints.size())
        *sup = *min_element(higherPoints.begin(),higherPoints.end());
    else *sup = t;
    if (lowerPoints.size())
        *inf = *max_element(lowerPoints.begin(),lowerPoints.end());
    else *inf = t;
    */
    // Or using trees
    m_trees[axis]->searchNode(t,inf,sup,false);
}

/******************************************************************************/


/************************* AI algo ********************************************/
// Comparer des points multivariés selon le temps d'attente dans la liste des voisins courants (nombre d'itération passé dans la liste m_curentNeighbours) 
bool alphaLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    return Utils::norm(nu->getAlpha(),2) < Utils::norm(mu->getAlpha(),2);
}
// Comparer des points multivariés selon le temps d'attente dans la liste des voisins courants (nombre d'itération passé dans la liste m_curentNeighbours) 
bool ageLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    return nu->getWaitingTime() < mu->getWaitingTime();
}
MultiVariatePointPtr<string> PiecewiseInterpolation::maxElement(int iteration)
{
    if (iteration%4)
        return *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),alphaLess);
    else
    {
        MultiVariatePointPtr<string> mu = *max_element(m_curentNeighbours.begin(), \
                                            m_curentNeighbours.end(),ageLess);
        for (MultiVariatePointPtr<string> nu : m_curentNeighbours)
            if (nu->alphaIsNull() && nu->getWaitingTime()==mu->getWaitingTime())
                return nu;
        return mu;
    }
}
MultiVariatePointPtr<string> PiecewiseInterpolation::getFirstMultivariatePoint()
{
    MultiVariatePointPtr<string> nu = make_shared<MultiVariatePoint<string>>(m_d,m_n,"");
    return nu;
}
void PiecewiseInterpolation::updateCurentNeighbours(MultiVariatePointPtr<string> nu)
{
 	// nu est le nouveau point d'interpolation (ce n'est plus un candidat)
  	m_curentNeighbours.remove(nu);

	// Incrémenter le temps d'attente de tous les voisins (candidats)
	for (MultiVariatePointPtr<string> mu : m_curentNeighbours)
  		mu->incrWaitingTime();

	// Ajouter les nouveaus candidats (voisins de nu)
	vector<string> childrenCodes;
  	for (int i=0; i<m_d; i++)
  	{
    	childrenCodes = BinaryTree::computeChildrenCodes((*nu)(i));
      	for (int k=0; k<int(childrenCodes.size()); k++)
      	{
        	MultiVariatePointPtr<string> mu = make_shared<MultiVariatePoint<string>>(*nu);
          	(*mu)(i) = childrenCodes[k];
          	if (isCorrectNeighbourToCurentPath(mu))
            	m_curentNeighbours.push_back(mu);
      	}
  	}
}
// Vérifier si index est déjà un point d'interpolation (ajouté dans path)
bool PiecewiseInterpolation::isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu)
{
    MultiVariatePoint<string> mu(*nu);
    string temp;
    for (int i=0; i<m_d; i++)
    {
        if ((*nu)(i).compare("")!=0)
        {
            temp = mu(i);
            mu(i) = BinaryTree::getParentCode(temp);
			// Si !indiceInPath(mu), nu ne peux pas être un candidat (il faut concerver un ensemble monotone)
            if (!indiceInPath(mu)) return false;
            mu(i) = temp;
        }
    }
    return true;
}
bool PiecewiseInterpolation::indiceInPath(MultiVariatePoint<string> index)
{
    for (MultiVariatePointPtr<string> nu : m_path)
    if (Utils::equals(index,*nu))
        return true;
    return false;
}
/******************************************************************************/

/*********************** Interpolation ****************************************/
// Fonction affine par morceaux
double PiecewiseInterpolation::basisFunction_1D(string code, double t, int axis)
{
    if (m_method==1)
    {
        if (code.compare("")==0) return 1;
        else if (code.compare("0")==0) return (t>0) ? 0 : -t;
        else if (code.compare("1")==0) return (t<0) ? 0 : t;
        else
        {
            double sup, inf, tcode = BinaryTree::getValueFromCode(code);
            computeBoundariesForBasisFunction(tcode,&inf,&sup,axis);
            if (t <= tcode && t >= inf) return ((t-inf)/(tcode-inf));
            else if (t >= tcode && t <= sup) return ((t-sup)/(tcode-sup));
            else return 0;
        }
    }
    else
    {
        if (code.compare("")==0) return 1;
        //else if (code.compare("0")==0) return (t>0) ? 0 : 0.5*t*(t - 1);
        //else if (code.compare("1")==0) return (t<0) ? 0 : 0.5*t*(t + 1);
        else if (code.compare("0")==0) return (t>0) ? 0 : -t*(t + 2);
        else if (code.compare("1")==0) return (t<0) ? 0 :  t*(2 - t);
        else
        {
            double sup, inf, tcode = BinaryTree::getValueFromCode(code);
            computeBoundariesForBasisFunction(tcode,&inf,&sup,axis);
            if (t <= sup && t >= inf) return (4/pow(sup-inf,2))*(sup-t)*(t-inf);
            else return 0;
        }
    }
}
/******************************************************************************/
