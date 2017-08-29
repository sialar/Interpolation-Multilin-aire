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
//
// 1: Concerne la méthode
//
// Core
// 1.0: Données d'entrée
		// Interpolé = fonction analytique implémenté dans la méthode evaluate de la classe Functions
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


// 1.1: Validation de la méthode
		// Taille de l'éspace des points de référence qui permet de valider la méthode
        int m_nbTestPoints;
        double m_infError;
        double m_mseError;
        vector<vector<double>> m_exactValues;
        vector<vector<double>> m_approxValues;
        vector<vector<double>> m_errors;

// Temps d'execution
		// Temps d'execution sans compter le cout de calcul de l'interpoléé en un point d'interpolation
        double m_runTime = 0.0;

// Sortie
		// Points d'interpolation sur chaque diréction
        vector<vector<double>> m_interpolationPoints;
		// Points multivariés d'interpolation
        vector<MultiVariatePoint<double>> m_interpolationNodes;
		// Points de test
        vector<MultiVariatePoint<double>> m_testPoints;

// Visualisation dans le terminal
        bool m_displayProgress = true;


    public:

        Interpolation() {};
        virtual ~Interpolation() {};

        Interpolation(FunctionsPtr f, int nIter);

		// Accesseurs
        const int nbEvals() { return m_nbEvals; };
        const double runTime() { return m_runTime; };
        const int maxIteration() { return m_maxIteration; };
        const int nbTestPoints() { return m_nbTestPoints; };
        void disableProgressDisplay() { m_displayProgress = false; };
        const vector<double> relativeErrors() { return m_infError; };
        const vector<MultiVariatePointPtr<T>>& path() { return m_path; };
        vector<MultiVariatePoint<double>> testPoints() { return m_testPoints; };
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        const vector<double>& interpolationPoints(int i) { return m_interpolationPoints[i]; };
        const vector<MultiVariatePoint<double>>& interpolationNodes() { return m_interpolationNodes; };
        const vector<vector<double>> parametersDomain() { return m_function->parametersDomain(); };
/******************************************************************************/
/************************ Points d'interpolation ******************************/
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

/******************************************************************************/
/*************************** Algorithme AI ************************************/
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


/******************************************************************************/
/************************** Fonctions de base *********************************/
        /**
          * Implémente la fonction de base (polynome de Lagrange global, fonction affine (ou quadratique) par morceaux (dépend de la version de l'algo AI)
          * \param code : ordre (indice ou code de Huffman) correspondant à une coordonnée d'un point d'interpolation
          * \param t : point d'évaluation (1d)
          * \param axis : direction concidéré (0 <= axis <= d-1)
          * \return la valeur au point t de la fonction de base correspondant au point 1d d'ordre code sur la direction axis
        */
        virtual double basisFunction_1D(T code, double t, int axis) = 0;

		// Construit la séquence de points de test aléatoirement
        void setRandomTestPoints(int nbTestPoints);
		// Lance l'algorithme d'interpolation adaptative
        void launchAIAlgo(bool debug);
		// Construit la séquence de points de test
        void setTestPoints(vector<MultiVariatePoint<double>> points);
		// Evaluer l'interpolé au point multivarié x
    vector<double> func(MultiVariatePoint<double> x);
    // Construire l'interpolé
    void setFunc(FunctionsPtr);
		// Calcule les résultats de l'interpolation
        void computeAIApproximationResults();
        void saveAIApproximationResults();
		// Calcule l'erreur d'interpolation au point d'interpolation courant (dernier choisi par l'algo AI)
        void computeLastAlphaNu(MultiVariatePointPtr<T> nu);
		// Construit l'interpolant en fonction des alpha et des fonctions de base
        vector<double> interpolation(MultiVariatePoint<double>& x, int end);
        vector<double> interpolation(MultiVariatePoint<double>& x);

        /************************* Other functions ****************************/

        void clearAll();
        void displayAll();
        void displayPath();
        void clearAllAlpha();
        void displayResults();
        void displayCurentNeighbours();
        void displayInterpolationPoints();
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
    auto start_time = chrono::steady_clock::now();

	// Initiatisation de m_curentNeighbours, m_path (séquence des points d'interpolation), argmax (le meilleur candidat)
    m_curentNeighbours.clear();
    m_path.clear();
    MultiVariatePointPtr<T> argmax = getFirstMultivariatePoint();
    MultiVariatePointPtr<T> old;
    m_curentNeighbours.push_back(argmax);
    int iteration = 0;

	// Si on atteint le nombre max d'itération, ou s'il n'y a plus de voisins possible
    while (!m_curentNeighbours.empty() && iteration < m_maxIteration)
    {
		// On commence par calculer (si ce n'est pas déjà fait) les alpha des candidats (points voisin à path)
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
		// On l'ajoute à m_path (ensemble courant des points d'interpolation)
        m_path.push_back(argmax);
        addInterpolationPoint(getPoint(argmax));

		// Mise à jour de la liste des candidats et du numéro de l'itération
        updateCurentNeighbours(argmax);
        iteration++;
    }

    auto end_time = chrono::steady_clock::now();
    std::chrono::duration<double> run_time = end_time - start_time;
    m_runTime = run_time.count();
}

template <typename T>
void Interpolation<T>::computeLastAlphaNu(MultiVariatePointPtr<T> nu)
{
    double basisFuncProd = 1.0;
    vector<double> res = func(getPoint(nu));
    m_nbEvals++;
    for (MultiVariatePointPtr<T> l : m_path)
    {
        basisFuncProd = 1.0;

        for (int p=0; p<m_d; p++)
            basisFuncProd *= basisFunction_1D((*l)(p),getPoint(nu)(p),p);
        for (int k=0; k<m_n; k++)
            res[k] -= l->getAlpha()[k] * basisFuncProd;
    }
    nu->setAlpha(res);
}

template <typename T>
vector<double> Interpolation<T>::interpolation(MultiVariatePoint<double>& x, int end)
{
  double l_prod = 1.0;
  vector<double> sum(m_n,0.0);
  for (int k=0; k<end; k++)
  {
      l_prod = 1.0;
      for (int i=0; i<m_d; i++)
          l_prod *= basisFunction_1D((*m_path[k])(i),x(i),i);
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
    m_function->displayParametersDomain();
	  cout << endl;
    displayInterpolationPoints();
	  cout << endl;
    displayResults();
}

template <typename T>
void Interpolation<T>::computeAIApproximationResults()
{
	// Implémenter ici le type d'erreur qu'on veut calculer
	vector<double> exactValue, approxValue, error;
	vector<double> exactValuesNorm, errorsNorm;

 	for (int i=0; i<m_nbTestPoints; i++)
	{
		exactValue = func(m_testPoints[i]);
		approxValue = interpolation(m_testPoints[i]);
		error = Utils::diff(exactValue,approxValue);

		m_exactValues.push_back(exactValue);
		m_approxValues.push_back(approxValue);
		m_errors.push_back(error);

		exactValuesNorm.push_back(Utils::norm(exactValue,2));
		errorsNorm.push_back(Utils::norm(error,2));
	}

	m_mseError = 0;
	vector<double> relativeErrors;
	double maxValue =  *max_element(exactValuesNorm.begin(),exactValuesNorm.end());

 	for (int i=0; i<m_nbTestPoints; i++)
	{
		relativeErrors.push_back(errorsNorm[i]*pow(10,5)/maxValue);
		m_mseError += pow(errorsNorm[i]*pow(10,5),2)/m_nbTestPoints;
	}

	m_infError = *max_element(relativeErrors.begin(), relativeErrors.end());
	m_mseError = sqrt(m_mseError);
}

template <typename T>
void Interpolation<T>::saveAIApproximationResults()
{
    cout << " - Approximation results are stored in data/output.dat" << endl;
    ofstream file("AI/data/output.dat", ios::out);
    if(file)
    {
        for (int i=0; i<int(m_testPoints.size()); i++)
        {
	          file << (i+1) << " ";
						for (int j=0; j<m_d; j++)
								file << /*setprecision(Utils::m_precision) <<*/ m_testPoints[i](j) << " ";
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

#endif