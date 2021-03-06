#include "../include/PiecewiseInterpolation.hpp"

PiecewiseInterpolation::PiecewiseInterpolation(FunctionsPtr f, int nIter, int method) :
    Interpolation(f, nIter), m_method(method)
{
    // Initialisation des arbres binaires sur chaque direction
    m_trees.resize(m_d);
    for (int i=0; i<m_d; i++)
      m_trees[i] = make_shared<BinaryTree>();
  }
void PiecewiseInterpolation::clearAllTrees()
{
    for (int i=0; i<int(m_trees.size()); i++)
        m_trees[i]->clearTree();
}

/******************************************************************************/
/************************ Points d'interpolation ******************************/
MultiVariatePoint<double> PiecewiseInterpolation::getPoint(MultiVariatePointPtr<string> nu)
{
    // Le point d'ordre nu est le point (x[nu(0)], .., x[nu(m_d)]), nu[i] est un code de Huffman
    MultiVariatePoint<double> point(m_d,0,0.0);
    for (int i=0; i<m_d; i++)
        // sur chaque direction, on calcule le point 1d correspondant au code de Huffman (*nu)(i)
        point(i) = BinaryTree::getValueFromCode((*nu)(i));
    return point;
}
void PiecewiseInterpolation::addInterpolationPoint(MultiVariatePoint<double>p)
{
    // Ajout dans m_interpolationNodes
    m_interpolationNodes.push_back(p);

    // Ajout de chaque coordonnée i dans m_interpolationPoints[i]
    // Ajout de chaque coordonnée i dans m_trees[i]
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

/******************************************************************************/
/*************************** Algorithme AI ************************************/
bool alphaLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    // Comparer des points multivariés selon le temps d'attente dans la liste des voisins courants (nombre d'itération passé dans la liste m_curentNeighbours)
    return Utils::norm(nu->getAlpha(),2) < Utils::norm(mu->getAlpha(),2);
}
bool ageLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    // Comparer des points multivariés selon le temps d'attente dans la liste des voisins courants (nombre d'itération passé dans la liste m_curentNeighbours)
    return nu->getNbWaitingIter() < mu->getNbWaitingIter();
}
MultiVariatePointPtr<string> PiecewiseInterpolation::maxElement(int iteration)
{
    // 1 fois sur 4, on choisi le point où l'erreur est la plus élevé
    if (iteration%4)
        return *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),alphaLess);
    else
    // 3 fois sur 4, on choisi le point qui a le plus attendu dans la liste des voisins
    // si il y en a plusieurs, on en choisit un où l'erreur est nulle pour éviter de bloquer une direction
    // En effet, si l'erreur en un point p vaut 0 elle le restaura si on ne procede que par adaptivité
    // Or parfois, une erreur peut être nulle par hasard (dépend de l'interpolé)
    // Exemple f(x,y) = sqrt(1-x**2)
    {
        MultiVariatePointPtr<string> mu = *max_element(m_curentNeighbours.begin(), \
                                            m_curentNeighbours.end(),ageLess);
        for (MultiVariatePointPtr<string> nu : m_curentNeighbours)
            if (nu->alphaIsNull() && nu->getNbWaitingIter()==mu->getNbWaitingIter())
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
  		mu->incrNbWaitingIter();

	// Ajouter les nouveaus candidats (voisins de nu)
	vector<string> childrenCodes;
  	for (int i=0; i<m_d; i++)
  	{
    	childrenCodes = BinaryTree::computeChildrenCodes((*nu)(i));
      	for (int k=0; k<int(childrenCodes.size()); k++)
      	{
        	 MultiVariatePointPtr<string> mu = make_shared<MultiVariatePoint<string>>(*nu);
           // mu est un voisin de nu si pour tout i mu(i) = nu(i) sauf pour un i0 où mu(i0) est un fils de nu(i0)
           (*mu)(i) = childrenCodes[k];
           // Vérifier que mu est un voisin de m_path courant
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
    // Vérifier si index est déjà un point d'interpolation (ajouté dans path)
    for (MultiVariatePointPtr<string> nu : m_path)
    if (Utils::equals(index,*nu))
        return true;
    return false;
}
/******************************************************************************/

/******************************************************************************/
/************************** Fonctions de base *********************************/
double PiecewiseInterpolation::basisFunction_1D(string code, double t, int axis)
{
    if (m_method==1)
    {
        // Fonction affine par morceaux (fonction chapeau)
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
        // Fonction quadratique par morceaux
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
void PiecewiseInterpolation::saveInterpolationBasisFunctions()
{
  ofstream file("AI/data/basis_functions.dat", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==1)
    {
      vector<double> x;
      for (int i=0; i<m_nbTestPoints; i++)
          x.push_back(m_testPoints[i](0));
      sort(x.begin(), x.end());
      file << m_path.size() << " " << m_nbTestPoints << " " << m_method << endl;
      MultiVariatePoint<double> p;
      for (int j=0; j<m_nbTestPoints; j++)
      {
        file << x[j];
        for (MultiVariatePointPtr<string> nu : m_path)
          file << " " <<  basisFunction_1D((*nu)(0),x[j],0);
        p = MultiVariatePoint<double>::toMonoVariatePoint(x[j]);
        file << " " <<  func(p)[0];
        file << endl;
      }
      for (MultiVariatePointPtr<string> nu : m_path)
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
void PiecewiseInterpolation::computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis)
{
    /*
    // Méthode 1:
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
    // Méthode 2: plus pérformante
    // Recherche des valeurs (inf et sup) les plus proches de t (de part et d'autre de t) dans l'arbre correspondant à la direction axis
    // Lors de la recherche, on ne concidère pas les nœuds descendants du nœud t. En d'autres termes on ne regarde que les nœuds dont la valeur correspondante (point 1d)
    // a un ordre plus petit que celui de t (point d'intepolation plus ancien)
    m_trees[axis]->searchNode(t,inf,sup);
}
