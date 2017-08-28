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

template <typename T>
class Interpolation
{
    public:
        static double m_precision;

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

// Temps d'execution
		// Temps d'execution sans compter le cout de calcul de l'interpoléé en un point d'interpolation
        double m_runTime = 0.0;
		// Temps total d'execution
		double m_totalTime = 0.0;
        chrono::time_point<chrono::_V2::steady_clock,chrono::duration<double>> m_lastCheckPt;

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
        const double totalTime() { return m_totalTime; };
        const int maxIteration() { return m_maxIteration; };
        const int nbTestPoints() { return m_nbTestPoints; };
        void disableProgressDisplay() { m_displayProgress = false; };
        const vector<double> relativeErrors() { return m_infError; };
        const vector<MultiVariatePointPtr<T>>& path() { return m_path; };
        vector<MultiVariatePoint<double>> testPoints() { return m_testPoints; };
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        const vector<double>& interpolationPoints(int i) { return m_interpolationPoints[i]; };
        const vector<MultiVariatePoint<double>>& interpolationNodes() { return m_interpolationNodes; };


		virtual void addInterpolationPoint(MultiVariatePoint<double> p) = 0;
		// Retourner le vecteur d'ordres du prochain point d'interpolation sélectionné par l'algo AI
		// 3 fois sur 4: on choisi celui ayant la plus grande erreur (algo AI)
		// 1 fois sur 4: on choisi celui qui a attendu le plus longtemps
		// Pour éviter de bloquer une direction dans le cas vaut l'erreur vaut par hasard 0
        virtual MultiVariatePointPtr<T> maxElement(int iteration) = 0;
		// Retourner le vecteur d'ordres du premier point d'interpolation (contient des codes de Huffman ou des indices)
        virtual MultiVariatePointPtr<T> getFirstMultivariatePoint() = 0;
		// Mettre à jour, en fonction du nouveau point d'interpolation, la liste des candidats courant
        virtual void updateCurentNeighbours(MultiVariatePointPtr<T> nu) = 0;
		// Vérifier si nu est bien un voisin de la séquence de points d'interpolation courante (vérifier la monotonie de l'ensemble)
        virtual bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<T> nu) = 0;
		// Implémente la fonction de base (polynome de Lagrange global, fonction affine (ou quadratique) par morceaux
        virtual double basisFunction_1D(T code, double t, int axis) = 0;
		// Méthode abstraite (le corp dépend de la vérsion concidérée (voir ligne 38))
		// Ajouter le point d'interpolation correspondant au vecteur d'ordre (multi-indice, ou vecteur de code de huffman selon la version de AI)
        virtual MultiVariatePoint<double> getPoint(MultiVariatePointPtr<T> nu) = 0;


		// Construit la séquence de points de test aléatoirement
        void setRandomTestPoints(int nbTestPoints);
		// Lance l'algorithme d'interpolation adaptative
        void launchAIAlgo(bool debug);
		// Construit la séquence de points de test
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };
		// Evaluer l'interpolé au point multivarié x
    vector<double> func(MultiVariatePoint<double> x);
    // Construire l'interpolé
    void setFunc(FunctionsPtr);
		// Calcule les résultats de l'interpolation
        void computeAIApproximationResults();
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
double Interpolation<T>::m_precision = numeric_limits<double>::digits10+1;

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
void Interpolation<T>::setRandomTestPoints(int nbTestPoints)
{
  m_nbTestPoints = nbTestPoints;
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
    m_lastCheckPt = chrono::steady_clock::now();

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
            cout << endl;
            Utils::separateur();
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
    std::chrono::duration<double> total_time = end_time - start_time;
    m_totalTime = total_time.count();
}

template <typename T>
void Interpolation<T>::computeLastAlphaNu(MultiVariatePointPtr<T> nu)
{
    double basisFuncProd = 1.0;
    auto stop_time = chrono::steady_clock::now();
    std::chrono::duration<double> delta = stop_time - m_lastCheckPt;
    m_runTime += delta.count();

    vector<double> res = func(getPoint(nu));
    m_lastCheckPt = chrono::steady_clock::now();
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
            cout << /*setprecision(m_precision) <<*/ *it << " ";
        cout << "}" << endl;
    }
}

template <typename T>
void Interpolation<T>::displayInterpolationMultiVariatePoints()
{
    cout << " - " << m_interpolationNodes.size() << "Interpolation nodes: { ";
    for (MultiVariatePoint<double> x : m_interpolationNodes)
        cout << setprecision(numeric_limits<double>::digits10+1) << x << " ";
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
        cout << " [" << *m_path[i] << ":" << setprecision(m_precision) << getPoint(m_path[i]);
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
        cout << "(" << (*nu) << ":" << setprecision(m_precision) << getPoint(nu) << ":";
        cout << Utils::vector2str(nu->getAlpha()) << ":" << nu << ") [" << nu->getWaitingTime() << "] | ";
    }
    cout << endl << endl;
}

template <typename T>
void Interpolation<T>::displayResults()
{
    cout << " - Interpolation error (pcm)" << endl;
    cout << " [e_inf] = " << m_infError << endl;
    cout << " [e_mse] = " << m_mseError << endl;
    cout << " - Number of calculation point = " << m_nbEvals << endl;
    cout << " - Total Time = " << m_totalTime << endl;
    cout << " - AI Run Time = " << m_runTime << endl;
}

template <typename T>
void Interpolation<T>::displayAll()
{
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
	vector<vector<double>> exactValues, approxValues, errors;
	vector<double> exactValuesNorm, errorsNorm;

 	for (int i=0; i<m_nbTestPoints; i++)
	{
		exactValue = func(m_testPoints[i]);
		approxValue = interpolation(m_testPoints[i]);
		error = Utils::diff(exactValue,approxValue);

		exactValues.push_back(exactValue);
		approxValues.push_back(approxValue);
		errors.push_back(error);

		exactValuesNorm.push_back(Utils::norm(exactValue,2));
		errorsNorm.push_back(Utils::norm(error,2));
	}

	m_mseError = 0;
	vector<double> relativeErrors;
	double maxValue =  Utils::max_elt(exactValuesNorm);
 	for (int i=0; i<m_nbTestPoints; i++)
	{
		relativeErrors.push_back(errorsNorm[i]*pow(10,5)/maxValue);
		m_mseError += pow(errorsNorm[i]*pow(10,5),2)/m_nbTestPoints;
	}

	m_infError = *max_element(relativeErrors.begin(), relativeErrors.end());
	m_mseError = sqrt(m_mseError);
}

#endif
