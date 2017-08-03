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
 *  \brief Classe générique, implémente l'algorithme d'Interpolation Adaptative
 *  \author SIALA Rafik
 *  \date 08/16
*/

template <typename T>
class Interpolation
{
    public:
        static double m_precision;

    protected:
        int m_d, m_n;
        int m_maxIteration;
        int m_nbTestPoints;
        int m_nbEvals = 0;
        int m_nbMethods = 3;

        double m_runTime = 0.0;
        double m_totalTime = 0.0;
        chrono::time_point<chrono::_V2::steady_clock,chrono::duration<double>> m_lastCheckPt;

        vector<vector<double>> m_realDomain;

        vector<vector<double>> m_interpolationPoints;
        vector<MultiVariatePoint<double>> m_interpolationNodes;
        vector<MultiVariatePoint<double>> m_testPoints;

        vector<MultiVariatePointPtr<T>> m_path;
        list<MultiVariatePointPtr<T>> m_curentNeighbours;

        map<string,map<method,vector<double>>> m_approxErrors;
        map<string,map<method,vector<double>>> m_approxResults;

        bool m_reactivity;
        bool m_displayProgress = true;
        FunctionsPtr m_function;

        vector<double> m_infErrors;
        vector<double> m_mseErrors;

    public:

        Interpolation() {};
        virtual ~Interpolation() {};
        /**
         *  Constructeur
         *  \param d: dimension de l'éspace de phase des paramètres
         *  \param core: Type de réacteur (MOX, UOX ou UOX-Gd)
         *  \param cs: Ensemble des sections efficaces à approcher
         *  \param nIter: Nombre d'itération dans l'algorithme AI (nombre de points d'interpolation)
        */
        Interpolation(int d, string core, vector<string> cs, int nIter);

        /************************* Data points ********************************/
        /**
         *  Retourner le nombre de points de calcul de l'interpolé
         *  \return Nombre d'évaluation
        */
        const int nbEvals() { return m_nbEvals; };
        /**
         *  Retourner le temps d'éxecution de l'algorithme AI (sans compter le cout des évaluations de l'interpolé)
         *  \return Temps d'exécution
        */
        const double runTime() { return m_runTime; };
        /**
         *  Retourner le temps d'éxecution de l'algorithme AI
         *  \return Temps d'exécution
        */
        const double totalTime() { return m_totalTime; };
        /**
         *  Retourner le nombre d'itérations dans l'algorithme AI
         *  \return Nombre d'itérations
        */
        const int maxIteration() { return m_maxIteration; };
        /**
         *  Retourner le nombre de points de référence (de tests)
         *  \return Nombre de points de référence
        */
        const int nbTestPoints() { return m_nbTestPoints; };
        /**
         *  Retourner les erreurs relatives absolues correspondants à chaque sections efficace
         *  \return Vecteur d'erreurs (l'élément i correspont à l'erreur commise lors de l'approximation de la section efficace d'indice i)
        */
        const vector<double> relativeErrors() { return m_infErrors; };
        /**
         *  Afficher la progression de la séquence de points d'interpolation et de la liste de priorité des voisins
        */
        void disableProgressDisplay() { m_displayProgress = false; };
        /**
        *  Retourner la séquence de multi-ordres des points d'interpolations
        *  \return Ensemble de multi-indices ou de vecteur de codes
        */
        const vector<MultiVariatePointPtr<T>>& path() { return m_path; };
        /**
        *  Methode virtuelle pour ajouter un nouveau point d'interpolation
        *  \param p Nouveau point choisi par AI
        */
        virtual void addInterpolationPoint(MultiVariatePoint<double> p) = 0;
        /**
        *  Retourner l'ensemble des points l'interpolation de chaque direction (chaque point est en dimension 1)
        *  \return Ensemble de points d'interpolation sur chaque variable
        */
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        /**
        *  Retourner le point correspondant au vecteur d'ordres nu
        *  \return Point d'interpolation ayant nu comme vecteur d'ordre
        */
        virtual MultiVariatePoint<double> getPoint(MultiVariatePointPtr<T> nu) = 0;
        /**
        *  Retourner l'ensemble des points d'interpolation multivariés
        *  \return Ensemble de points d'interpolation
        */
        const vector<MultiVariatePoint<double>>& interpolationNodes() { return m_interpolationNodes; };
        /**
        *  Retourner l'ensemble des points l'interpolation sur la direction i
        *  \return Ensemble de points discrétisation dans la variable i
        */
        const vector<double>& interpolationPoints(int i) { return m_interpolationPoints[i]; };

        /**
        *  Créer un vecteur de points multivariés aléatoirement
        *  \param nbTestPoints Nombre de points multivariés
        */
        void setRandomTestPoints(int nbTestPoints);
        /**
        *  Retourner l'ensemble de points de référence
        *  \return Ensemble de points de test
        */
        vector<MultiVariatePoint<double>> testPoints() { return m_testPoints; };
        /**
        *  Créer la grille de référence à partir vecteur entré en paramètre
        *  \param points Ensemble de points de test
        */
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };

        /**
        *  Définir le type de réacteur
        *  \param c Type de réacteur (MOX, UOX ou UOX-Gd)
        */
        void setFunc(string c);
        /**
        *  Définir le type de réacteur et les sections efficaces
        *  \param c Type de réacteur (MOX, UOX ou UOX-Gd)
        *  \param c Type de sections efficaces (macro_totale0, macro_fission0, macro_absorption0, ...)
        */
        void setFunc(string c, vector<string> cs);
        /**
        *  Retourner la valeur de l'interpolé au point x
        *  \return x Point d'évaluation
        */
        vector<double> func(MultiVariatePoint<double> x);


        /************************* AI algo ************************************/
        /**
        *  Calculer la réactivité obtenu avec l'algorithme AI
        */
        void computeReactivity();
        /**
        *  Calculer les facteur de multiplication efficace K_eff (en chaque point de test) en utilsant la méthode m
        *  \param m Méthode d'approximation utilisée (Cocagne, Tucker ou AI)
        *  \return les facteurs de multiplication efficace en chaque point de référence
        */
        vector<double> computeKinf(method m);
        /**
        *  Calculer les résultats d'approximations ainsi que les erreurs obtenues avec la méthode parcimonieuse.
        *  Puis les stocker dans les attributs m_approxResults et m_approxErrors
        */
        void computeAIApproximationResults();
        /**
        *  Calculer l'erreur d'interpolation au point multivarié courant (condidat ayant nu comme veteur d'ordres) dans l'exécution de l'algorithme AI.
        *  \param nu Vecteur d'ordres du condidat courant
        */
        void computeLastAlphaNu(MultiVariatePointPtr<T> nu);
        /**
        *  Méthode qui implémente l'algorithme d'Interpolation Adaptative.
        *  \param debug Boolean pour éxecuter avec (ou sans) affichage de la progression de l'algo
        */
        void launchAIAlgo(bool debug);
        /**
        *  Méthode implémentée dan les classes filles (dépend du type d'ordre).
        *  Elle permet de choisir un point parmi les condidat courant (les points voisin à la séquence courante de points d'interpolation).
        *  Le choix est fait d'une manière adaptative dans la plupart des fois.
        *  Parfois (25% des cas), on choisit le condidat qui attendu le plus pour être sur de bien explorer tout l'éspace.
        *  \return Le nouveau point d'interpolation
        */
        virtual MultiVariatePointPtr<T> maxElement(int iteration) = 0;
        /**
        *  Créer le vecteur d'ordres du premier point d'interpolation ( [0,0,0,0,0] ou ["","","","",""]).
        *  \return Le premier point d'interpolation
        */
        virtual MultiVariatePointPtr<T> getFirstMultivariatePoint() = 0;
        /**
        *  Mettre à jour la liste des condidats, une fois le point d'ordre nu ajouté à la séquence d'interpolation.
        *  \return Dernier point d'interpolation ajouté
        */
        virtual void updateCurentNeighbours(MultiVariatePointPtr<T> nu) = 0;
        /**
        *  Vérifier si le nouveau condidat est bien voisin à la séquence courante de points d'interpolation.
        *  Cette vérification permet de s'assurer de la monotonie de la séquence
        *  \param nu Nouveau condidat potentiel
        *  \return true si nu correspont à un condidat valide, false sinon
        */
        virtual bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<T> nu) = 0;

        /************************* Interpolation ******************************/
        /**
        *  Calculer la valuer de la fonction de base (polynome de Lagrange ou fonction quadratique par morceaux) associée à "code" au point "t" sur la direction "axis"
        *  \param code Ordre (1d) du point d'interpolation correspondant à la fonction de base en question
        *  \param t Point d'évaluation (1d)
        *  \param axis Direction
        *  \return Valeur de la fonction de base
        */
        virtual double basisFunction_1D(T code, double t, int axis) = 0;
        /**
        *  Calculer l'approximation des sections efficaces au point de référence x en utilsant "end" itérations
        *  \param x Point de référence
        *  \param end Nombre d'itérations utilisé
        *  \return Approximation des sections efficaces en x
        */
        vector<double> interpolation(MultiVariatePoint<double>& x, int end);
        /**
        *  Calculer l'approximation des sections efficaces au point de référence x en utilsant le nombre maximal d'itérations
        *  \param x Point de référence
        *  \return Approximation des sections efficaces en x
        */
        vector<double> interpolation(MultiVariatePoint<double>& x);

        /************************* Other functions ****************************/
        void comp();
        void clearAll();
        void displayAll();
        void displayPath();
        void saveResults();
        void clearAllAlpha();
        void displayResults();
        void displayRealDomain();
        void readDataAndResults();
        void saveReactivityInFile();
        void readReactivityFromFile();
        void saveTuckerResultsInFile();
        void displayCurentNeighbours();
        void displayCrossSectionNames();
        void displayInterpolationPoints();
        void readReferencePointsFromFile();
        void saveApproximationResultsInFile();
        void readApproximationResultsFromFile();
        void computeTuckerApproximationResults();
        void saveInterpolationMultiVariatePoints();
        void displayInterpolationMultiVariatePoints();
};

template <typename T>
using InterpolationPtr = shared_ptr<Interpolation<T>>;

template <typename T>
double Interpolation<T>::m_precision = numeric_limits<double>::digits10+1;

template <typename T>
Interpolation<T>::Interpolation(int d, string core, vector<string> cs, int nIter)
{
    m_d = d;
    m_n = cs.size();
    m_realDomain.resize(d);
    m_maxIteration = nIter;
    m_interpolationPoints.resize(d);
    m_function = make_shared<Functions>(core, cs);
    m_reactivity = m_function->reactivityIsComputable();
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
void Interpolation<T>::setFunc(string c, vector<string> vr)
{
    m_n = vr.size();
    m_function->setCoreType(c);
    m_function->setCrossSectionType(vr);
}

template <typename T>
void Interpolation<T>::setFunc(string c)
{
    m_function->setCoreType(c);
    m_function->setAllCrossSectionType();
    m_n = m_function->getCrossSections().size();
}

template <typename T>
vector<double> Interpolation<T>::func(MultiVariatePoint<double> x)
{
    vector<double> fast_res, res;
    for (int i=0; i<m_d; i++)
        x(i) = Utils::convertToFunctionDomain(m_realDomain[i][0], m_realDomain[i][1], x(i));
    fast_res = m_function->fast_evaluate(x);
    return fast_res;
}

/*************************** Data Points **************************************/
template <typename T>
void Interpolation<T>::setRandomTestPoints(int nbTestPoints)
{
  vector<MultiVariatePoint<double>> testPoints;
  testPoints.resize(nbTestPoints);
  for (int j=0; j<nbTestPoints; j++)
      testPoints[j] = Utils::createRandomMultiVariatePoint(m_d);
  setTestPoints(testPoints);
}
/******************************************************************************/

/***************************** AI algo ****************************************/
template <typename T>
void Interpolation<T>::launchAIAlgo(bool debug)
{
    auto start_time = chrono::steady_clock::now();
    m_lastCheckPt = chrono::steady_clock::now();

    m_curentNeighbours.clear();
    m_path.clear();
    MultiVariatePointPtr<T> argmax = getFirstMultivariatePoint();
    MultiVariatePointPtr<T> old;
    m_curentNeighbours.push_back(argmax);
    int iteration = 0;

    while (!m_curentNeighbours.empty() && iteration < m_maxIteration)
    {
        for (MultiVariatePointPtr<T> nu : m_curentNeighbours)
            if (!nu->alphaAlreadyComputed())
                computeLastAlphaNu(nu);

        if (debug)
        {
            cout << endl;
            Utils::separateur();
            displayPath();
            displayCurentNeighbours();
        }

        argmax = maxElement(iteration);
        m_path.push_back(argmax);
        addInterpolationPoint(getPoint(argmax));
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
/******************************************************************************/

/*************************** Interpolation ************************************/
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
/******************************************************************************/

/*************************** Display function *********************************/
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
void Interpolation<T>::displayRealDomain()
{
    cout << " - Real domain of the cross section function : ";
    for (int i=0; i<m_d-1; i++)
      cout << "[" << m_realDomain[i][0] << "," << m_realDomain[i][1] << "] x ";
    cout << "[" << m_realDomain[m_d-1][0] << "," << m_realDomain[m_d-1][1] << "]" << endl;
}

template <typename T>
void Interpolation<T>::displayCrossSectionNames()
{
    cout << " - Cross section types : [ " << m_function->getCrossSections()[0];
    for (int i=1; i<m_n; i++)
        cout  << " ; " << m_function->getCrossSections()[i];
    cout << " ]" << endl;
}

template <typename T>
void Interpolation<T>::displayResults()
{
    vector<double> co_err_inf, tu_err_inf, ai_err_inf;
    vector<double> co_err_mse, tu_err_mse, ai_err_mse;
    for (int i=0; i<m_n; i++)
    {
        string csName = m_function->getCrossSections()[i];
        co_err_inf.push_back(Utils::maxAbsValue(m_approxErrors[csName][Cocagne]));
        tu_err_inf.push_back(Utils::maxAbsValue(m_approxErrors[csName][Tucker]));
        ai_err_inf.push_back(Utils::maxAbsValue(m_approxErrors[csName][AI]));
        co_err_mse.push_back(Utils::computeMseError(m_approxResults[csName][Apollo],m_approxResults[csName][Cocagne]));
        tu_err_mse.push_back(Utils::computeMseError(m_approxResults[csName][Apollo],m_approxResults[csName][Tucker]));
        ai_err_mse.push_back(Utils::computeMseError(m_approxResults[csName][Apollo],m_approxResults[csName][AI]));
    }
    cout << " - Interpolation error using Cocagne (pcm)" << endl;
    cout << " [e_inf] = ";
    Utils::displayValues(co_err_inf);
    cout << " [e_mse] = ";
    Utils::displayValues(co_err_mse);
    cout << " - Interpolation error using Tucker (pcm)" << endl;
    cout << " [e_inf] = ";
    Utils::displayValues(tu_err_inf);
    cout << " [e_mse] = ";
    Utils::displayValues(tu_err_mse);
    cout << " - Interpolation error using AI (pcm)" << endl;
    cout << " [e_inf] = ";
    Utils::displayValues(ai_err_inf);
    cout << " [e_mse] = ";
    Utils::displayValues(ai_err_mse);

    if (m_reactivity)
    {
      cout << " - Reactivity error using Cocagne (pcm) = " << Utils::maxAbsValue(m_approxErrors["reactivity"][Cocagne]) << " " << Utils::computeMseError(m_approxResults["reactivity"][Apollo],m_approxResults["reactivity"][Cocagne]) << endl;
      cout << " - Reactivity error using Tucker (pcm) = " << Utils::maxAbsValue(m_approxErrors["reactivity"][Tucker]) << " " << Utils::computeMseError(m_approxResults["reactivity"][Apollo],m_approxResults["reactivity"][Tucker]) << endl;
      cout << " - Reactivity error using AI (pcm) = " << Utils::maxAbsValue(m_approxErrors["reactivity"][AI]) << " " << Utils::computeMseError(m_approxResults["reactivity"][Tucker],m_approxResults["reactivity"][AI]) << endl;
    }

    cout << " - Number of calculation point = " << m_nbEvals << endl;
    cout << " - Total Time = " << m_totalTime << endl;
    cout << " - AI Run Time = " << m_runTime << endl;
}

template <typename T>
void Interpolation<T>::displayAll()
{
    displayRealDomain();
    displayCrossSectionNames();
    Utils::separateur();
    displayInterpolationPoints();
    Utils::separateur();
    displayResults();
}
/******************************************************************************/

/*************************** Results ******************************************/
template <typename T>
void Interpolation<T>::computeAIApproximationResults()
{
    for (int i=0; i<m_n; i++)
    {
        string csName = m_function->getCrossSections()[i];
        vector<double> ai_res, ai_err;
        double val, realValue, maxValue = Utils::maxAbsValue(m_approxResults[csName][Apollo]);
        for (int j=0; j<m_nbTestPoints; j++)
        {
            val = interpolation(m_testPoints[j])[i];
            realValue = m_approxResults[csName][Apollo][j];
            ai_res.push_back(val);
            ai_err.push_back(pow(10,5)*(val-realValue)/maxValue);
        }
        m_approxResults[csName].insert(pair<method,vector<double>>(AI,ai_res));
        m_approxErrors[csName].insert(pair<method,vector<double>>(AI,ai_err));

        m_infErrors.push_back(Utils::maxAbsValue(m_approxErrors[csName][AI]));
        m_mseErrors.push_back(Utils::computeMseError(m_approxResults[csName][Apollo],m_approxResults[csName][AI]));
    }
}

template <typename T>
void Interpolation<T>::computeTuckerApproximationResults()
{
    for (int i=0; i<m_n; i++)
    {
        string csName = m_function->getCrossSections()[i];
        vector<double>  tu_err, tu_res;
        double val, realValue, maxValue = Utils::maxAbsValue(m_approxResults[csName][Apollo]);
        for (int j=0; j<m_nbTestPoints; j++)
        {
            MultiVariatePoint<double> x(m_testPoints[j]);
            for (int k=0; k<m_d; k++)
                x(k) = Utils::convertToFunctionDomain(m_realDomain[k][0], m_realDomain[k][1], m_testPoints[j](k));
            val = m_function->fast_evaluate(x)[i];
            realValue = m_approxResults[csName][Apollo][j];
            tu_res.push_back(val);
            tu_err.push_back(pow(10,5)*(val-realValue)/maxValue);
        }
        m_approxResults[csName].insert(pair<method,vector<double>>(Tucker_bis,tu_res));
        m_approxErrors[csName].insert(pair<method,vector<double>>(Tucker_bis,tu_err));
    }
}

template <typename T>
vector<double> Interpolation<T>::computeKinf(method m)
{
    double nufi0, nufi1, tot0, tot1, scatt000, scatt001, scatt010, scatt011;
    double num, denom;
    vector<double> k_eff;
    for (int i=0; i<m_nbTestPoints; i++)
    {
        nufi0 = m_approxResults["macro_nu*fission0"][m][i];
        nufi1 = m_approxResults["macro_nu*fission1"][m][i];
        tot0 = m_approxResults["macro_totale0"][m][i];
        tot1 = m_approxResults["macro_totale1"][m][i];
        scatt000 = m_approxResults["macro_scattering000"][m][i];
        scatt001 = m_approxResults["macro_scattering001"][m][i];
        scatt010 = m_approxResults["macro_scattering010"][m][i];
        scatt011 = m_approxResults["macro_scattering011"][m][i];
        num = nufi0 * (tot1 - scatt011) + nufi1*scatt010;
        denom = (tot0 - scatt000 ) * (tot1 - scatt011) - scatt001*scatt010;
        k_eff.push_back(num / denom);
    }
    return k_eff;
}

template <typename T>
void Interpolation<T>::computeReactivity()
{
    vector<double> kinfApollo = computeKinf(Apollo);
    vector<double> kinfAI = computeKinf(AI);
    vector<double> reactivityError(m_nbTestPoints);
    for (int i=0; i<m_nbTestPoints; i++)
        reactivityError[i] = (1/kinfApollo[i] - 1/kinfAI[i]) * pow(10,5);
    m_approxResults["reactivity"].insert(pair<method,vector<double>>(AI,kinfAI));
    m_approxErrors["reactivity"].insert(pair<method,vector<double>>(AI,reactivityError));
}

template <typename T>
void Interpolation<T>::readDataAndResults()
{
    readReferencePointsFromFile();
    readApproximationResultsFromFile();
    if (m_reactivity) readReactivityFromFile();
}

template <typename T>
void Interpolation<T>::saveResults()
{
    computeAIApproximationResults();
    saveApproximationResultsInFile();
    if (m_reactivity)
    {
        computeReactivity();
        saveReactivityInFile();
    }
    saveInterpolationMultiVariatePoints();
}

template <typename T>
void Interpolation<T>::readReferencePointsFromFile()
{
    string csName = m_function->getCrossSections()[0];
    string s = Utils::replace(csName,"*","_");
    ifstream file(m_function->realDataDirPath() + "/FinalResults/" + s, ios::in);
    if(file)
    {
        string line;
        vector<double> data;
        vector<double> max(m_d,-numeric_limits<double>::max());
        vector<double> min(m_d,numeric_limits<double>::max());
        MultiVariatePoint<double> p(5,0,0);
        while (getline(file, line))
        {
            data = Utils::str2vector(line);
            for (int i=0; i<5; i++)
                p(i) = data[i+1];
            m_testPoints.push_back(p);
            for (int i=0; i<m_d; i++)
            {
                if (p(i) > max[i]) max[i] = p(i);
                if (p(i) < min[i]) min[i] = p(i);
            }
        }
        for (int i=0; i<m_d; i++)
        {
            m_realDomain[i].push_back(min[i]);
            m_realDomain[i].push_back(max[i]);
        }
        m_nbTestPoints = m_testPoints.size();
        for (int i=0; i<m_nbTestPoints; i++)
            for (int j=0; j<m_d; j++)
                m_testPoints[i](j) = Utils::adaptCoordsToFunctionDomain(min[j], max[j], m_testPoints[i](j));
        file.close();
    }
    else cerr << "Error while opening " << csName << " file!" << endl;
}

template <typename T>
void Interpolation<T>::readApproximationResultsFromFile()
{
    for (int k=0; k<m_n; k++)
    {
        string csName = m_function->getCrossSections()[k];
        string s = Utils::replace(csName,"*","_");
        ifstream file(m_function->realDataDirPath() + "/FinalResults/" + s, ios::in);
        if(file)
        {
            string line;
            vector<double> data, ap_res, co_res, tu_res;
            vector<double> co_err, tu_err;
            while (getline(file, line))
            {
                data = Utils::str2vector(line);
                ap_res.push_back(data[6]);
                co_res.push_back(data[7]);
                tu_res.push_back(data[8]);
                co_err.push_back(data[9]);
                tu_err.push_back(data[10]);
            }
            map<method,vector<double>> mapRes, mapErr;
            mapRes.insert(pair<method,vector<double>>(Apollo,ap_res));
            mapRes.insert(pair<method,vector<double>>(Cocagne,co_res));
            mapRes.insert(pair<method,vector<double>>(Tucker,tu_res));
            m_approxResults.insert(pair<string,map<method,vector<double>>>(csName,mapRes));
            mapErr.insert(pair<method,vector<double>>(Cocagne,co_err));
            mapErr.insert(pair<method,vector<double>>(Tucker,tu_err));
            m_approxErrors.insert(pair<string,map<method,vector<double>>>(csName,mapErr));
            file.close();
        }
        else cerr << "Error while opening " << csName << " file!" << endl;
    }

}

template <typename T>
void Interpolation<T>::saveApproximationResultsInFile()
{
    for (string csName : m_function->getCrossSections())
    {
        string s = Utils::replace(csName,"*","_");
        ofstream file(Utils::projectPath + "AI/data/" + m_function->getCoreType() + "/" + s, ios::out);
        if(file)
        {
            for (int i=0; i<m_nbTestPoints; i++)
            {
              MultiVariatePoint<double> x(m_testPoints[i]);
              for (int j=0; j<m_d; j++)
                  x(j) = Utils::convertToFunctionDomain(m_realDomain[j][0], m_realDomain[j][1], m_testPoints[i](j));

              file << (i+1) << " ";
              for (int j=0; j<m_d; j++)
                  file << setprecision(m_precision) << x(j) << " ";
              file << setprecision(m_precision) << m_approxResults[csName][Apollo][i] << " ";
              file << setprecision(m_precision) << m_approxResults[csName][Cocagne][i] << " ";
              file << setprecision(m_precision) << m_approxResults[csName][Tucker][i] << " ";
              file << setprecision(m_precision) << m_approxResults[csName][AI][i] << " ";
              file << setprecision(m_precision) << m_approxErrors[csName][Cocagne][i] << " ";
              file << setprecision(m_precision) << m_approxErrors[csName][Tucker][i] << " ";
              file << setprecision(m_precision) << m_approxErrors[csName][AI][i] << endl;

            }
            file.close();
        }
        else cerr << "Error while opening the file!" << endl;
    }
}

template <typename T>
void Interpolation<T>::saveTuckerResultsInFile()
{
    computeTuckerApproximationResults();
    for (int k=0; k<m_n; k++)
    {
        string csName = m_function->getCrossSections()[k];
        string s = Utils::replace(csName,"*","_");
        ofstream file(Utils::projectPath + "AI/data/RegeneratedResults/" + m_function->getCoreType() + "/" + s, ios::out);
        if(file)
        {
            for (int i=0; i<m_nbTestPoints; i++)
            {
              MultiVariatePoint<double> x(m_testPoints[i]);
              for (int j=0; j<m_d; j++)
                  x(j) = Utils::convertToFunctionDomain(m_realDomain[j][0], m_realDomain[j][1], m_testPoints[i](j));

              //file << (i+1) << " ";
              file << setprecision(m_precision) << m_approxResults[csName][Tucker][i] << " ";
              file << setprecision(m_precision) << m_approxResults[csName][Tucker_bis][i] << endl;
            }
            file.close();
        }
        else cerr << "Error while opening the file!" << endl;
    }
}

template <typename T>
void Interpolation<T>::readReactivityFromFile()
{
    ifstream file(m_function->realDataDirPath() + "/FinalResults/ReactivityError", ios::in);
    if(file)
    {
        string line;
        vector<double> data, ap_rea, co_rea, tu_rea;
        vector<double> co_err, tu_err;
        while (getline(file, line))
        {
            data = Utils::str2vector(line);
            ap_rea.push_back(data[6]);
            co_rea.push_back(data[7]);
            tu_rea.push_back(data[8]);
            co_err.push_back(data[9]);
            tu_err.push_back(data[10]);
        }
        map<method,vector<double>> mapRes, mapErr;
        mapRes.insert(pair<method,vector<double>>(Apollo,ap_rea));
        mapRes.insert(pair<method,vector<double>>(Cocagne,co_rea));
        mapRes.insert(pair<method,vector<double>>(Tucker,tu_rea));
        m_approxResults.insert(pair<string,map<method,vector<double>>>("reactivity",mapRes));
        mapErr.insert(pair<method,vector<double>>(Cocagne,co_err));
        mapErr.insert(pair<method,vector<double>>(Tucker,tu_err));
        m_approxErrors.insert(pair<string,map<method,vector<double>>>("reactivity",mapErr));
        file.close();
    }
    else cerr << "Error while opening the file!" << endl;
}

template <typename T>
void Interpolation<T>::saveReactivityInFile()
{
    ofstream file(Utils::projectPath + "AI/data/" + m_function->getCoreType() + "/ReactivityError", ios::out );
    if(file)
    {
        for (int i=0; i<m_nbTestPoints; i++)
        {
          MultiVariatePoint<double> x(m_testPoints[i]);
          for (int j=0; j<m_d; j++)
              x(j) = Utils::convertToFunctionDomain(m_realDomain[j][0], m_realDomain[j][1], m_testPoints[i](j));

          file << (i+1) << " ";
          for (int j=0; j<m_d; j++)
              file << setprecision(m_precision) << x(j) << " ";

          file << setprecision(m_precision) << m_approxResults["reactivity"][Apollo][i] << " ";
          file << setprecision(m_precision) << m_approxResults["reactivity"][Cocagne][i] << " ";
          file << setprecision(m_precision) << m_approxResults["reactivity"][Tucker][i] << " ";
          file << setprecision(m_precision) << m_approxResults["reactivity"][AI][i] << " ";
          file << setprecision(m_precision) << m_approxErrors["reactivity"][Cocagne][i] << " ";
          file << setprecision(m_precision) << m_approxErrors["reactivity"][Tucker][i] << " ";
          file << setprecision(m_precision) << -m_approxErrors["reactivity"][AI][i] << endl;
        }
        file.close();
    }
    else cerr << "Error while opening the file!" << endl;
}
/******************************************************************************/

template <typename T>
void Interpolation<T>::comp()
{
    for (string csName : m_function->getCrossSections())
    {
        vector<double> data, indices, saved_res;
        string s = Utils::replace(csName,"*","_");

        ifstream file_in(Utils::projectPath + "AI/data/" + m_function->getCoreType() + "/" + s, ios::out);
        if(file_in)
        {
            string line;
            while (getline(file_in, line))
            {
                data = Utils::str2vector(line);
                indices.push_back(data[0]);
                saved_res.push_back(data[8]);
            }
            file_in.close();
        }
        else cerr << "Error while opening the file!" << endl;

        ofstream file(Utils::projectPath + "AI/data/" + m_function->getCoreType() + "/comp_" + s, ios::out);
        if(file)
        {
            for (int i=0; i<m_nbTestPoints; i++)
            {
                MultiVariatePoint<double> x(m_testPoints[i]);
                for (int j=0; j<m_d; j++)
                    x(j) = Utils::convertToFunctionDomain(m_realDomain[j][0], m_realDomain[j][1], m_testPoints[i](j));

                file << indices[i] << " ";
                file << setprecision(m_precision) << saved_res[i] << " ";
                file << setprecision(m_precision) << m_function->fast_evaluate(x)[0] << " ";
                file << setprecision(m_precision) << m_function->evaluate(x)[0] << endl;
            }
            file.close();
        }
        else cerr << "Error while opening the file!" << endl;
    }
}

/******************************************************************************/

template <typename T>
void Interpolation<T>::saveInterpolationMultiVariatePoints()
{

  ofstream file(Utils::projectPath + "AI/data/interpolation_points.dat", ios::out);
  if(file)
  {
      double t;
      file << m_interpolationNodes.size() << endl;
      for (MultiVariatePoint<double> x : m_interpolationNodes)
      {
          for (int i=0; i<m_d; i++)
          {
              t = Utils::convertToFunctionDomain(m_realDomain[i][0], m_realDomain[i][1], x(i));
              file << setprecision(numeric_limits<double>::digits10+1) << t << " ";
          }
          file << endl;
      }
      file.close();
  }
  else cerr << "Error while opening the file!" << endl;

}

#endif
