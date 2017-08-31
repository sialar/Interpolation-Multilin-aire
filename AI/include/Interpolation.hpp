#ifndef INTERPOLATION
#define INTERPOLATION

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <memory>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <limits>

#include "MultiVariatePoint.hpp"
#include "BinaryTree.hpp"
#include "Functions.hpp"
#include "Utils.hpp"

using namespace std;

/**
 *  \file Interpolation.hpp
 *  \brief Classe générique abstraite qui implémente l'algorithme d'interpolation adaptative
 *  \author SIALA Rafik
 *  \date 08/16
*/
template <typename T>
class Interpolation
{
    protected:
/******************************************************************************/
/**************************** Données d'entrée ********************************/
      // Fonction qu'on veut interpoler, elle est soit de type:
      // AnalyticalFunctions: dans ce cas on veut valider la méthode en approchant une fonctions analytique implmenté dans la méthode evaluate
      // RealDataFunctions: dans ce cas on veut tester la méthode sur des données réelles
      FunctionsPtr m_function;
      // m_d : dimension de l'éspace de départ de l'interpolé
		  // m_n : dimension de l'éspace d'arrivé de l'interpolé
		  // exemple: f: (x_1, .., x_{m_d}) --> (y_1, .., y_{m_n})
      int m_d, m_n;
		  // Nombre d'itérations dans l'algo AI (= Nombre de points d'interpolation)
      int m_maxIteration;
		  // Nombre total de calculs de f efféctués (Nombre d'appel à evaluate de Functions: en pratique, correspond au nombre d'appels à un code coûteux)
		  // Attention: m_nbEvals != m_maxIteration, car m_nbEvals comprend les points qui seront évalués mais qui ne seront points des points d'interpolation
      int m_nbEvals = 0;
		  // Nombre de versions de l'algo AI (Points de Leja + Polynomes de Lagrange), (Arbre binaire + Fonctions affines par morceaux) ou (Arbre binaire + Fonctions quad par morceaux)
      int m_nbMethods = 3;
		  // Chemin correspondant à la séquence ordonnée de points d'interpolation, construit par l'algorithme AI
      vector<MultiVariatePointPtr<T>> m_path;
		  // Ensemble de voisins (candidats) dans l'itération courante de l'algorithme AI
      list<MultiVariatePointPtr<T>> m_curentNeighbours;
      // 0 si LagrangeInterpolation, 1 ou 2 si PiecewiseInterpolation

/******************************************************************************/
/************************ Validation de la méthode ****************************/
		    // Taille de l'éspace des points de référence qui permet de valider la méthode
        int m_nbTestPoints;
        // Erreur relative absolue calculée sur l'ensemble de points de référence
        double m_infError;
        // Erreur quadratique moyenne calculée sur l'ensemble de points de référence
        double m_mseError;
        // Ensemble des valeurs exacte de f sur la liste de points m_testPoints
        vector<vector<double>> m_exactValues;
        // Ensemble des valeurs approchés de f sur la liste de points m_testPoints
        vector<vector<double>> m_approxValues;
        // Ensemble des erreurs d'interpolation calculées sur la liste de points m_testPoints
        vector<vector<double>> m_errors;

/******************************************************************************/
/**************************** Sortie de l'algo ********************************/
        // Temps d'execution de la'lgorithme AI (calculé dans la méthode launchAIAlgo())
        double m_runTime = 0.0;
    		//0Ensemble des points d'interpolation sur chaque diréction
        vector<vector<double>> m_interpolationPoints;
		    // Ensemble des points multivariés d'interpolation
        vector<MultiVariatePoint<double>> m_interpolationNodes;
		    // Ensemble de points de test (points de référence) sur lequel on évalue la méthode
        // Dans le cas de fonctions analytiques, il s'agit d'une liste aléatoire de points
        // Dans le cas de données réelles, il s'agit de la liste des points lue à partir du fichier de données exemple "croos_section.dat"
        vector<MultiVariatePoint<double>> m_testPoints;

    public:

        Interpolation() {};
        virtual ~Interpolation() {};
        /**
          * Constructeur
          * \param f : fonction à interpoler
          * \param nIter : nombre d'itérations dans l'algorithme AI
        */
        Interpolation(FunctionsPtr f, int nIter);
        /**
          * Supprimer toutes les données
        */
        void clearAll();
        /**
          * Sauvegarder tous les résultats de l'interpolation
        */
        void saveAll();

/******************************************************************************/
/**************************** Données d'entrée ********************************/
        /**
          * Accesseur à l'attribut m_nbEvals
          * \return m_nbEvals, nombre de points de calcul
        */
        const int nbEvals() { return m_nbEvals; };
        /**
          * Accesseur à l'attribut m_runTime
          * \return m_runTime, temps d'éxecution de l'algo AI
        */
        const double runTime() { return m_runTime; };
        /**
          * Accesseur à l'attribut m_maxIteration
          * \return m_maxIteration, nombre d'itérations dans l'algo AI
        */
        const int maxIteration() { return m_maxIteration; };
        /**
          * Accesseur à l'attribut m_nbTestPoints
          * \return m_nbTestPoints, nombre de points de référence (points de test)
        */
        const int nbTestPoints() { return m_nbTestPoints; };
        /**
          * Evaluer l'interpolé au point multivarié x (appelle la méthode evaluate de Functions)
          * \param x : point multivarié
          * \return valeur de l'interpolé f en x
        */
        vector<double> func(MultiVariatePoint<double> x);
        /**
          * Mutateur de l'attribut m_function
          * \param f : fonction qu'on veut interpoler
        */
        void setFunc(FunctionsPtr f);
        /**
          * Retourne le domaine de définition de l'interpolé f
          * \return domaine de définition de f (attribut m_parametersDomain de m_function)
        */
        const vector<vector<double>> parametersDomain() { return m_function->parametersDomain(); };

/******************************************************************************/
/************************ Validation de la méthode ****************************/
        /**
          * Accesseur à l'attribut m_infError
          * \return m_infError, valeur de l'erreur relative absolue
        */
        const vector<double> relativeErrors() { return m_infError; };
        /**
          * Accesseur à l'attribut m_testPoints
          * \return m_testPoints, ensemble de points de test pour évaluer la méthode
        */
        vector<MultiVariatePoint<double>> testPoints() { return m_testPoints; };
        /**
          * Construit la séquence de points de test, m_testPoints, aléatoirement
          * \param nbTestPoints : nombre de points de test
        */
        void setRandomTestPoints(int nbTestPoints);
        /**
          * Mutateur de l'attribut m_testPoints
          * \param points : ensemble de points multivariés de test
        */
        void setTestPoints(vector<MultiVariatePoint<double>> points);
        /**
          * Calculer les résultats d'approximation obtenus par l'algo AI
        */
        void computeAIApproximationResults();
        /**
          * Sauvegarder les résultats d'approximation obtenus par l'algo AI dans le fichier data/output.dat
        */
        void saveAIApproximationResults();

/******************************************************************************/
/************************ Points d'interpolation ******************************/
        /**
          * Retourne l'ensemble de points de discretisation sur la direction i
          * \param i : direction
          * \return vecteur de points de discretisation sur la direction i
        */
        const vector<double>& interpolationPoints(int i) { return m_interpolationPoints[i]; };
        /**
          * Accesseur à l'attribut m_interpolationNodes
          * \return m_interpolationNodes, nsemble des points multivariés d'interpolation
        */
        const vector<MultiVariatePoint<double>>& interpolationNodes() { return m_interpolationNodes; };
        /**
          * Retourne le point multivarié correspondant au vecteur d'ordres nu
          * \param nu : vecteur d'ordres indices ou codes de Huffman (dépend du type générique T)
          * \return point multivarié correspondant à nu
        */
        virtual MultiVariatePoint<double> getPoint(MultiVariatePointPtr<T> nu) = 0;
        /**
          * Ajoute le nouveau point d'interpolation p
          * \param p : nouveau point d'interpolation
        */
		    virtual void addInterpolationPoint(MultiVariatePoint<double> p) = 0;
        /**
        * Stocker les valeurs des fonctions de bases dans basis_functions.dat
        */
        virtual void saveInterpolationBasisFunctions() = 0;
        /**
          * Stocker les points d'interpolation dans le fichier data/interpolation_points.dat
        */
        void saveInterpolationPoints();


/******************************************************************************/
/*************************** Algorithme AI ************************************/
        /**
          * Accesseur à l'attribut m_path
          * \return m_path, chemin correspondant à la séquence ordonnée de points d'interpolation, construit par l'algo AI
        */
        const vector<MultiVariatePointPtr<T>>& path() { return m_path; };
        /**
        * Retourner le premier vecteur d'ordres dans l'algo AI
        * \return le premier vecteur d'ordres
        */
        virtual MultiVariatePointPtr<T> getFirstMultivariatePoint() = 0;
        /**
          * Retourner le vecteur d'ordres du prochain point d'interpolation sélectionné par l'algo AI
          * 3 fois sur 4: on choisi celui ayant la plus grande erreur (algo AI)
      		* 1 fois sur 4: on choisi celui qui a attendu le plus longtemps pour éviter de bloquer une direction dans le cas où l'erreur vaut 0
          * \param iteration : numéro de l'itération courante
          * \return teration : vecteur d'ordres du prochain point d'interpolation
        */
        virtual MultiVariatePointPtr<T> maxElement(int iteration) = 0;
        /**
          * Vérifier si le vecteur d'ordre index existe dans m_path (correspond à un point d'interpolation)
          * \param index : vecteur d'ordre
          * \return true si et seulement si index est dans m_path
        */
        virtual bool indiceInPath(MultiVariatePoint<T> index) = 0;
        /**
          * Mettre à jour, en fonction du nouveau point d'interpolation correspondant au vecteur d'ordre nu, la liste des candidats courant
          * \param nu : nouveau vecteur d'ordre choisi par l'algo AI
        */
        virtual void updateCurentNeighbours(MultiVariatePointPtr<T> nu) = 0;
        /**
          * Vérifier si nu est bien un voisin de la séquence de points d'interpolation courante (vérifier la monotonie de l'ensemble de points choisis)
          * \param nu : vecteur d'ordres
          * \return true si et seulement si nu est voisin à l'ensemble de points de m_path
        */
        virtual bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<T> nu) = 0;
        /**
          * Implémente l'algorithme d'interpolation adaptative
          * \param debug : si true, affiche sur le terminal la progression de m_path et m_currentNeighbors
        */
        void launchAIAlgo(bool debug);
        /**
          * Réinitialiser tous les alpha
        */
        void clearAllAlpha();
        /**
          * Calculer l'erreur d'interpolation alpha au point d'interpolation courant
          * \param nu : vecteur d'ordre qui correspond au dernier point d'interpolation choisi par l'algo AI
        */
        void computeLastAlphaNu(MultiVariatePointPtr<T> nu);
        /**
          * Calculer la valeur appproché de f au point multivarié x en utilisant end iterations dans l'algo AI
          * \param x : point multivarié
          * \param end : nombre d'itération à utiliser pour la construction de l'interpolant
          * \return valeur vectorielle approchée de f obtenue par l'algo AI
        */
        vector<double> interpolation(MultiVariatePoint<double>& x, int end);
        /**
          * Calculer la valeur appproché de f, au point multivarié x, obtenue au bout de m_maxIteration iterations
          * \param x : point multivarié
          * \return valeur vectorielle approchée de f obtenue par l'algo AI
        */
        vector<double> interpolation(MultiVariatePoint<double>& x);
        /**
        * Implémente la fonction de base (polynome de Lagrange global, fonction affine (ou quadratique) par morceaux (dépend de la version de l'algo AI)
        * \param code : ordre (indice ou code de Huffman) correspondant à une coordonnée d'un point d'interpolation
        * \param t : point d'évaluation (1d)
        * \param axis : direction concidéré (0 <= axis <= d-1)
        * \return la valeur au point t de la fonction de base correspondant au point 1d d'ordre code sur la direction axis
        */
        virtual double basisFunction_1D(T code, double t, int axis) = 0;
        /**
          * Stocker dans data/interpolation_progression.dat la valeur de l'approximation à chaque iteration sur tous les points de test
        */
        void saveInterpolationProgression();
        /**
          * Stocker les élement de m_path (points d'interpolations) dans l'ordre dans le fichier data/path.dat
        */
        void savePath();


/******************************************************************************/
/************************ Fonctions d'affichage *******************************/
        /**
          * Afficher tous les résultats (domaine de définition de f, l'ensemble de points construit par l'algo AI, les erreurs d'interpolation ainsi que
          * le temps d'execution) obtenus par l'algo AI
        */
        void displayAll();
        /**
          * Afficher l'ensemble courant de points d'interpolation dans l'ordre
        */
        void displayPath();
        /**
          * Afficher les erreurs d'interpolation, le nombre points de calcul ainsi que le temps d'execution
        */
        void displayResults();
        /**
          * Afficher l'ensemble de voisins courant à m_path
        */
        void displayCurentNeighbours();
        /**
          * Afficher l'ensemble des points de test direction par direction
        */
        void displayInterpolationPoints();
        /**
          * Afficher l'ensemble des points multivariés de test
        */
        void displayInterpolationMultiVariatePoints();
};

template <typename T>
using InterpolationPtr = shared_ptr<Interpolation<T>>;

template <typename T>
Interpolation<T>::Interpolation(FunctionsPtr f, int nIter)
{
    m_d = f->getD();
    m_n = f->getN();
    m_function = f;
    m_maxIteration = nIter;
    m_interpolationPoints.resize(m_d);
}

template <typename T>
void Interpolation<T>::clearAll()
{
    m_path.clear();
    m_curentNeighbours.clear();
    m_interpolationNodes.clear();
    for (int i=0; i<m_d; i++)
        m_interpolationPoints[i].clear();
}


template <typename T>
void Interpolation<T>::clearAllAlpha()
{
    for (MultiVariatePointPtr<T> nu : m_path)
        nu->reinit();
}

template <typename T>
void Interpolation<T>::saveAll()
{
    savePath();
    saveInterpolationPoints();
    saveAIApproximationResults();
    saveInterpolationProgression();
    saveInterpolationBasisFunctions();
}

template <typename T>
vector<double> Interpolation<T>::func(MultiVariatePoint<double> x)
{
    return m_function->evaluate(x);
}

template <typename T>
void Interpolation<T>::setTestPoints(vector<MultiVariatePoint<double>> points)
{
  m_testPoints = points;
  m_nbTestPoints = points.size();
}

template <typename T>
void Interpolation<T>::setRandomTestPoints(int nbTestPoints)
{
  vector<MultiVariatePoint<double>> testPoints;
  testPoints.resize(nbTestPoints);
  for (int j=0; j<nbTestPoints; j++)
      testPoints[j] = Utils::createRandomMultiVariatePoint(m_d);
  setTestPoints(testPoints);
}

template <typename T>
void Interpolation<T>::launchAIAlgo(bool debug)
{
    // Début du calcul du temps d'execution
    auto start_time = chrono::steady_clock::now();

	  // Initiatisation de m_curentNeighbours, m_path (séquence des points d'interpolation), argmax (le meilleur candidat courant)
    m_curentNeighbours.clear();
    m_path.clear();
    // Constuit le vecteur d'ordres correspondant au premier point d'interpolation
    MultiVariatePointPtr<T> argmax = getFirstMultivariatePoint();
    // Ajouter le premier point à la liste des voisins courant
    m_curentNeighbours.push_back(argmax);
    // Initilisation de la valeur de l'iteration courante
    int iteration = 0;

	  // On s'arrete si on atteint le nombre max d'itération, ou s'il n'y a plus de voisins possible
    while (!m_curentNeighbours.empty() && iteration < m_maxIteration)
    {
		    // On commence par calculer (si ce n'est pas déjà fait) les alpha des candidats (points voisins à path)
        for (MultiVariatePointPtr<T> nu : m_curentNeighbours)
            if (!nu->alphaAlreadyComputed())
                computeLastAlphaNu(nu);

		    // Pour suivre la progression de l'algo sur le terminal
        if (debug)
        {
            cout << endl << endl;
            displayPath();
            displayCurentNeighbours();
        }

		    // Une fois tous les alpha calculés, on choisi le meilleur
        argmax = maxElement(iteration);
		    // On ajoute le vecteur d'ordres correspondant dans m_path (ensemble courant des points d'interpolation)
        m_path.push_back(argmax);
        // On ajoute le point d'interpolation correspondant dans m_interpolationNodes et m_interpolationPoints
        addInterpolationPoint(getPoint(argmax));

		    // Mise à jour de la liste des candidats
        updateCurentNeighbours(argmax);
        // Mise à jour de la valeur de l'iteration courante
        iteration++;
    }

    // Fin du calcul du temps d'execution
    auto end_time = chrono::steady_clock::now();
    std::chrono::duration<double> run_time = end_time - start_time;
    m_runTime = run_time.count();
}

template <typename T>
void Interpolation<T>::computeLastAlphaNu(MultiVariatePointPtr<T> nu)
{
    // alpha(nu) = f(y_{nu}) - I_{Lambda_{nu}}(f)(y_{nu})
    // ou encore: alpha = f(y_{nu}) - \sum_{mu<nu} alpha(mu) * H_{mu}(y_{nu}) (rapport page 12)
    double basisFuncProd = 1.0;
    vector<double> res = func(getPoint(nu)); // f(y_{nu})
    // Incrémentation du nombre de calcul de f éffectué par l'algo AI
    m_nbEvals++;
    // représente la somme sur les points d'interpolation courants
    for (MultiVariatePointPtr<T> mu : m_path)
    {
        // Calcul de H, produit tensoriel (H_{mu}(f)(y_{nu})) des fonctions hiérarchiques au point d'interpolation y_{nu}
        basisFuncProd = 1.0;
        for (int p=0; p<m_d; p++)
            basisFuncProd *= basisFunction_1D((*mu)(p),getPoint(nu)(p),p);

        for (int k=0; k<m_n; k++)
            // correspond à alpha(mu) * H_{mu}(f)(y_{nu})
            res[k] -= mu->getAlpha()[k] * basisFuncProd;
    }
    // Sauvegarder la valeur de l'erreur au point d'ordre nu
    nu->setAlpha(res);
}

template <typename T>
vector<double> Interpolation<T>::interpolation(MultiVariatePoint<double>& x, int end)
{
  // I_{Lambda_{n}}(f)(x) = \sum_{k<n} alpha(nu_k) * H_{nu_k}(x) (rapport page 12)
  double l_prod = 1.0;
  vector<double> sum(m_n,0.0);
  for (int k=0; k<end; k++)
  {
      // Calcul de H_{nu_k}, produit tensoriel des fonctions hiérarchiques au point x
      l_prod = 1.0;
      for (int i=0; i<m_d; i++)
          l_prod *= basisFunction_1D((*m_path[k])(i),x(i),i);
      // Calucl de la somme des H_{nu_k} multipliés par les alpha des k premiers points
      for (int i=0; i<m_n; i++)
          sum[i] += (m_path[k]->getAlpha())[i] * l_prod;
  }
  return sum;
}

template <typename T>
vector<double> Interpolation<T>::interpolation(MultiVariatePoint<double>& x)
{
  return interpolation(x,m_maxIteration);
}

template <typename T>
void Interpolation<T>::displayInterpolationPoints()
{
    vector<double>::iterator it;
    for (int i=0; i<m_d; ++i)
    {
        cout << " - " << m_interpolationPoints[i].size() << " points in direction " << i << " : { ";
        for (it=m_interpolationPoints[i].begin(); it!=m_interpolationPoints[i].end(); it++)
            cout << /*setprecision(Utils::m_precision) <<*/ *it << " ";
        cout << "}" << endl;
    }
}

template <typename T>
void Interpolation<T>::displayInterpolationMultiVariatePoints()
{
    cout << " - " << m_interpolationNodes.size() << " interpolation nodes: { ";
    for (MultiVariatePoint<double> x : m_interpolationNodes)
        cout << /*setprecision(Utils::m_precision) <<*/ x << " ";
    cout << "}" << endl;
}

template <typename T>
void Interpolation<T>::displayPath()
{
    // Format: (nu:points[nu]) --> () ...
    cout << " - Path =";
    int n = m_path.size();
    for (int i=0; i<n; i++)
    {
        if (i>0) cout << "\t ";
        cout << " " << i << " :";
        cout << " [" << *m_path[i] << ":" << setprecision(Utils::m_precision) << getPoint(m_path[i]);
        cout << ":" << m_path[i] << "] --> alpha" << *m_path[i] << " = " << Utils::vector2str(m_path[i]->getAlpha()) << endl;
    }
    cout << endl;
}

template <typename T>
void Interpolation<T>::displayCurentNeighbours()
{
    cout << " - Curent neighbours (" << m_curentNeighbours.size() << ") = ";
    for (MultiVariatePointPtr<T> nu : m_curentNeighbours)
    {
        cout << "(" << (*nu) << ":" << setprecision(Utils::m_precision) << getPoint(nu) << ":";
        cout << Utils::vector2str(nu->getAlpha()) << ":" << nu << ") [" << nu->getNbWaitingIter() << "] | ";
    }
    cout << endl << endl;
}

template <typename T>
void Interpolation<T>::displayResults()
{
    cout << " - Interpolation error (pcm):";
    cout << " [e_inf] = " << m_infError << " | [e_mse] = " << m_mseError << endl;
    cout << endl << " - Number of calculation point = " << m_nbEvals << endl;
    cout << endl << " - AI run time = " << m_runTime << endl;
}

template <typename T>
void Interpolation<T>::displayAll()
{
    cout << endl;
    // Afficher le domaine de définition de f
    m_function->displayParametersDomain();
	  cout << endl;
    // Afficher l'ensemble final de points d'interpolation construit par l'algo AI
    displayInterpolationPoints();
	  cout << endl;
    // Afficher les erreurs d'interpolation ainsi que le temps de calcul
    displayResults();
}

template <typename T>
void Interpolation<T>::computeAIApproximationResults()
{
	// Calculer les erreurs d'interpolation (erreur relative absolue, erreur quadratique moyenne) sur
  // la grille des points de test, m_testPoints
	vector<double> exactValue, approxValue, error;
	vector<double> exactValuesNorm, errorsNorm;

 	for (int i=0; i<m_nbTestPoints; i++)
	{
    // valeur exacte de f au ième point de référence
		exactValue = func(m_testPoints[i]);
    // valeur approché de f au ième point de référence
		approxValue = interpolation(m_testPoints[i]);
    // erreur d'interpolation de f au ième point de référence
		error = Utils::diff(exactValue,approxValue);

		m_exactValues.push_back(exactValue);
		m_approxValues.push_back(approxValue);
		m_errors.push_back(error);

    // Calcul des normes des valeurs exacte pour la recherche du sup d'entre eux
		exactValuesNorm.push_back(Utils::norm(exactValue,2));
    // Calcul des normes des erreurs
		errorsNorm.push_back(Utils::norm(error,2));
	}

	m_mseError = 0;
	vector<double> relativeErrors;
  // max des valeurs exacte, pour le calcul de l'erreur relative
	double maxValue =  *max_element(exactValuesNorm.begin(),exactValuesNorm.end());

 	for (int i=0; i<m_nbTestPoints; i++)
	{
    // calcul et sauvegarde des erreurs relative sur chaque point de test
		relativeErrors.push_back(errorsNorm[i]*pow(10,5)/maxValue);
    // calcul de l'erreur quadratique moyenne
		m_mseError += pow(errorsNorm[i]*pow(10,5),2)/m_nbTestPoints;
	}

	m_infError = *max_element(relativeErrors.begin(), relativeErrors.end());
	m_mseError = sqrt(m_mseError);
}

template <typename T>
void Interpolation<T>::saveAIApproximationResults()
{
    cout << " - Approximation results are stored in AI/data/output.dat" << endl;
    // Sauvegarde des résultats de l'interpolation de f
    ofstream file("AI/data/output.dat", ios::out);
    if(file)
    {
        for (int i=0; i<int(m_testPoints.size()); i++)
        {
            MultiVariatePoint<double> x(m_testPoints[i]);
            for (int j=0; j<m_d; j++)
                x(j) = Utils::convertToFunctionDomain(parametersDomain()[j][0], parametersDomain()[j][1], m_testPoints[i](j));

	          file << (i+1) << " ";
						for (int j=0; j<m_d; j++)
								file << setprecision(Utils::m_precision) << x(j) << " ";
            file << "| ";
						for (int j=0; j<m_n; j++)
								file << /*setprecision(Utils::m_precision) <<*/ m_exactValues[i][j] << " ";
            file << "| ";
            for (int j=0; j<m_n; j++)
                file << /*setprecision(Utils::m_precision) <<*/ m_approxValues[i][j] << " ";
            file << "| ";
            for (int j=0; j<m_n; j++)
                file << /*setprecision(Utils::m_precision) <<*/ m_errors[i][j] << endl;
				}
        file.close();
    }
        else cerr << "Error while opening the file!" << endl;
}

template <typename T>
void Interpolation<T>::saveInterpolationPoints()
{

  ofstream file("AI/data/interpolation_points.dat", ios::out);
  if(file)
  {
      double t;
      file << m_interpolationNodes.size() << endl;
      for (MultiVariatePoint<double> x : m_interpolationNodes)
      {
          for (int i=0; i<m_d; i++)
          {
              t = Utils::convertToFunctionDomain(parametersDomain()[i][0], parametersDomain()[i][1], x(i));
              file << setprecision(numeric_limits<double>::digits10+1) << t << " ";
          }
          file << endl;
      }
      file.close();
  }
  else cerr << "Error while opening the file!" << endl;
}

template <typename T>
void Interpolation<T>::saveInterpolationProgression()
{
  ofstream file("AI/data/interpolation_progression.dat", ios::out | ios::trunc);
  if(file)
  {
      if (m_d==1)
      {
          vector<double> x;
          for (int i=0; i<m_nbTestPoints; i++)
              x.push_back(m_testPoints[i](0));
          sort(x.begin(), x.end());
          MultiVariatePoint<double> p;
          vector<double> tempPath;
          for (int j=0; j<m_nbTestPoints; j++)
          {
              for (int i=0; i<int(m_path.size()); i++)
              {
                  p = MultiVariatePoint<double>::toMonoVariatePoint(x[j]);
                  file << interpolation(p, i+1)[0] << " ";
              }
              file << endl;
          }
      }
      file.close();
  }
  else
      cerr << "Error while opening the file!" << endl;
}

template <typename T>
void Interpolation<T>::savePath()
{
  ofstream file("AI/data/path.dat", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==2 || m_d==3)
    {
        for (int i=0; i<m_d; i++)
            file << m_interpolationPoints[i].size() << " ";
        file << endl << m_path.size() << endl;
        MultiVariatePoint<double> x(m_d,0,0.0);
        for (MultiVariatePointPtr<T> nu : m_path)
        {
            x = getPoint(nu);
            for (int i=0; i<m_d; i++)
                file << x(i) << " ";
            file << endl;
        }
        for (int i=0; i<int(m_path.size())-1; i++)
            file << Utils::vector2str(m_path[i]->getAlpha()) << " ; ";
        file << Utils::vector2str(m_path[m_path.size()-1]->getAlpha()) << endl;
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}
#endif
