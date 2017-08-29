#include "../include/MixedInterpolation.hpp"

MixedInterpolation::MixedInterpolation(FunctionsPtr f, int nIter, MultiVariatePoint<int> methods) :
    Interpolation(f, nIter)
{
    // Charger la séquence de points de Leja 1d depuis le fichier data/leja_sequence.dat
    m_lejaSequence = Utils::loadLejaSequenceFromFile(m_maxIteration);
    m_trees.resize(m_d);
    m_methods = methods;
    // Initialisation des arbres binaires sur chaque direction
    for (int i=0; i<m_d; i++)
      m_trees[i] = make_shared<BinaryTree>();
}
void MixedInterpolation::clearAllTrees()
{
    for (int i=0; i<int(m_trees.size()); i++)
        m_trees[i]->clearTree();
}

/******************************************************************************/
/************************ Points d'interpolation ******************************/
MultiVariatePoint<double> MixedInterpolation::getPoint(MultiVariatePointPtr<string> nu)
{
    // Le point d'ordre nu est le point (x[nu(0)], .., x[nu(m_d)]), nu[i] est soit un code de Huffman soit un entier correspondant à un indice dans la séquence de Leja
    MultiVariatePoint<double> point(m_d,0,0.0);
    for (int i=0; i<m_d; i++)
    {
        if (m_methods(i))
            // sur chaque direction (de version 1 ou 2), on calcule le point 1d correspondant au code de Huffman (*nu)(i)
            point(i) = BinaryTree::getValueFromCode((*nu)(i));
        else
            // sur chaque direction i, on prend le point ième point 1d de la séquence de Leja
            point(i) = m_lejaSequence[stoi((*nu)(i))];
    }
    return point;
}
void MixedInterpolation::addInterpolationPoint(MultiVariatePoint<double>p)
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
            if (m_methods(i)) m_trees[i]->addNode(p(i));
        }
    }
}

/******************************************************************************/
/*************************** Algorithme AI ************************************/
bool customAlphaLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    // Comparer des points multivariés selon le temps d'attente dans la liste des voisins courants (nombre d'itération passé dans la liste m_curentNeighbours)
    return Utils::norm(nu->getAlpha(),2) < Utils::norm(mu->getAlpha(),2);
}
bool customAgeLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    // Comparer des points multivariés selon le temps d'attente dans la liste des voisins courants (nombre d'itération passé dans la liste m_curentNeighbours)
    return nu->getNbWaitingIter() < mu->getNbWaitingIter();
}
MultiVariatePointPtr<string> MixedInterpolation::maxElement(int iteration)
{
    // 1 fois sur 4, on choisi le point où l'erreur est la plus élevé
    if (iteration%4)
        return *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),customAlphaLess);
    else
    // 3 fois sur 4, on choisi le point qui a le plus attendu dans la liste des voisins
    // si il y en a plusieurs, on en choisit un où l'erreur est nulle pour éviter de bloquer une direction
    // En effet, si l'erreur en un point p vaut 0 elle le restaura si on ne procede que par adaptivité
    // Or parfois, une erreur peut être nulle par hasard (dépend de l'interpolé)
    // Exemple f(x,y) = sqrt(1-x**2)
    {
        MultiVariatePointPtr<string> mu = *max_element(m_curentNeighbours.begin(), \
                                            m_curentNeighbours.end(),customAgeLess);
        for (MultiVariatePointPtr<string> nu : m_curentNeighbours)
            if (nu->alphaIsNull() && nu->getNbWaitingIter()==mu->getNbWaitingIter())
                return nu;
        return mu;
    }
}
MultiVariatePointPtr<string> MixedInterpolation::getFirstMultivariatePoint()
{
    MultiVariatePointPtr<string> nu = make_shared<MultiVariatePoint<string>>(m_d,m_n,"");
    for (int i=0; i<m_d; i++)
        if (!m_methods(i))
            (*nu)(i) = "0";
    return nu;
}
void MixedInterpolation::updateCurentNeighbours(MultiVariatePointPtr<string> nu)
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
        if (m_methods(i))
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
        else
        {
            MultiVariatePointPtr<string> mu = make_shared<MultiVariatePoint<string>>(*nu);
            // mu est un voisin de nu si pour tout i mu(i) = nu(i) sauf pour un i0 où nu1(i0) = nu(i0)+1
            (*mu)(i) = to_string(stoi((*mu)(i))+1);
            // Vérifier que mu est un voisin de m_path courant
            if (isCorrectNeighbourToCurentPath(mu))
                m_curentNeighbours.push_back(mu);
        }
    }
}
bool MixedInterpolation::isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu)
{
    MultiVariatePoint<string> mu(*nu);
    string temp;
    for (int i=0; i<m_d; i++)
    {
        if (m_methods(i))
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
        else
        {
            if ((*nu)(i).compare("0")!=0)
            {
                mu(i) = to_string(stoi(mu(i))-1);
				        // Si !indiceInPath(mu), nu ne peux pas être un candidat (il faut concerver un ensemble monotone)
	            if (!indiceInPath(mu)) return false;
                mu(i) = to_string(stoi(mu(i))+1);
            }
        }
    }
    return true;
}
bool MixedInterpolation::indiceInPath(MultiVariatePoint<string> index)
{
    // Vérifier si index est déjà un point d'interpolation (ajouté dans path)
    for (MultiVariatePointPtr<string> nu : m_path)
    if (Utils::equals(index,*nu))
        return true;
    return false;
}

/******************************************************************************/
/************************** Fonctions de base *********************************/
double MixedInterpolation::basisFunction_1D(string code, double t, int axis)
{
    if (m_methods(axis)==0)
    {
        // Polynôme hiérarchique de Lagrange
        if (code.compare("0")==0) return 1;
        int k = stoi(code);
        double prod = 1;
        for (int i=0; i<k; ++i)
            prod *= (t-m_lejaSequence[i]) / (m_lejaSequence[k]-m_lejaSequence[i]) ;
        return prod;
    }
    else if (m_methods(axis)==1)
    {
        // Fonnctions affine par morceaux (fonction chapeau)
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
        // Fonnctions quadratique par morceaux
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
void MixedInterpolation::computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis)
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
/******************************************************************************/
