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
        int m_d;
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

        map<method,vector<double>> m_approxErrors;
        map<method,vector<double>> m_approxResults;

        bool m_displayProgress = true;
        FunctionsPtr m_function;

    public:

        Interpolation() {};
        virtual ~Interpolation() {};
        Interpolation(int d, string core, string cs, int nIter);

        /************************* Data points ********************************/
        const int nbEvals() { return m_nbEvals; };
        const double runTime() { return m_runTime; };
        const double totalTime() { return m_totalTime; };
        const int maxIteration() { return m_maxIteration; };
        const int nbTestPoints() { return m_nbTestPoints; };
        void disableProgressDisplay() { m_displayProgress = false; };
        const vector<MultiVariatePointPtr<T>>& path() { return m_path; };
        virtual void addInterpolationPoint(MultiVariatePoint<double> p) = 0;
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        virtual MultiVariatePoint<double> getPoint(MultiVariatePointPtr<T> nu) = 0;
        const vector<MultiVariatePoint<double>>& interpolationNodes() { return m_interpolationNodes; };
        const vector<double>& interpolationPoints(int i) { return m_interpolationPoints[i]; };

        void setRandomTestPoints(int nbTestPoints);
        vector<MultiVariatePoint<double>> testPoints() { return m_testPoints; };
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };

        void setFunc(string c, string cs);
        FunctionsPtr getFunc() { return m_function; };
        double func(MultiVariatePoint<double> x);


        /************************* AI algo ************************************/
        void computeReactivity();
        double computeKinf(method m);
        void computeAIApproximationResults();
        void computeLastAlphaNu(MultiVariatePointPtr<T> nu);
        int launchAIAlgo(bool debug);
        virtual MultiVariatePointPtr<T> maxElement(int iteration) = 0;
        virtual MultiVariatePointPtr<T> getFirstMultivariatePoint() = 0;
        virtual void updateCurentNeighbours(MultiVariatePointPtr<T> nu) = 0;
        virtual bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<T> nu) = 0;

        /************************* Interpolation ******************************/
        virtual double basisFunction_1D(T code, double t, int axis) = 0;
        double interpolation(MultiVariatePoint<double>& x, int end);
        double interpolation(MultiVariatePoint<double>& x);

        /************************* Other functions ****************************/
        void clearAll();
        void displayAll();
        void displayPath();
        void saveResults();
        void clearAllAlpha();
        void displayResults();
        void displayFunction();
        void displayRealDomain();
        void readDataAndResults();
        void displayCurentNeighbours();
        void displayInterpolationPoints();
        void readReferencePointsFromFile();
        void saveApproximationResultsInFile();
        void readApproximationResultsFromFile();
};

template <typename T>
using InterpolationPtr = shared_ptr<Interpolation<T>>;

template <typename T>
double Interpolation<T>::m_precision = numeric_limits<double>::digits10+1;

template <typename T>
Interpolation<T>::Interpolation(int d, string core, string cs, int nIter)
{
    m_d = d;
    m_realDomain.resize(d);
    m_maxIteration = nIter;
    m_interpolationPoints.resize(d);
    m_function = make_shared<Functions>(core, cs);
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
void Interpolation<T>::setFunc(string c, string cs)
{
    m_function->setCoreType(c);
    m_function->setCrossSectionType(cs);
}

template <typename T>
double Interpolation<T>::func(MultiVariatePoint<double> x)
{
    for (int i=0; i<m_d; i++)
        x(i) = Utils::convertToFunctionDomain(m_realDomain[i][0], m_realDomain[i][1], x(i));
    return m_function->evaluate(x);
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
int Interpolation<T>::launchAIAlgo(bool debug)
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
    return iteration;
}

template <typename T>
void Interpolation<T>::computeLastAlphaNu(MultiVariatePointPtr<T> nu)
{
    double basisFuncProd = 1.0;
    auto stop_time = chrono::steady_clock::now();
    std::chrono::duration<double> delta = stop_time - m_lastCheckPt;
    m_runTime += delta.count();

    double res = func(getPoint(nu));
    m_lastCheckPt = chrono::steady_clock::now();
    m_nbEvals++;
    for (MultiVariatePointPtr<T> l : m_path)
    {
        basisFuncProd = 1.0;
        for (int p=0; p<m_d; p++)
            basisFuncProd *= basisFunction_1D((*l)(p),getPoint(nu)(p),p);
        res -= l->getAlpha() * basisFuncProd;
    }
    nu->setAlpha(res);
}
/******************************************************************************/

/*************************** Interpolation ************************************/
template <typename T>
double Interpolation<T>::interpolation(MultiVariatePoint<double>& x, int end)
{
  double l_prod = 1.0, sum = 0.0;
  for (int k=0; k<end; k++)
  {
      l_prod = 1.0;
      for (int i=0; i<m_d; i++)
          l_prod *= basisFunction_1D((*m_path[k])(i),x(i),i);
      sum += m_path[k]->getAlpha() * l_prod;
  }
  return sum;
}

template <typename T>
double Interpolation<T>::interpolation(MultiVariatePoint<double>& x)
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
        cout << ":" << m_path[i] << "] --> alpha" << *m_path[i] << " = " << m_path[i]->getAlpha() << endl;
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
        cout << nu->getAlpha() << ":" << nu << ") [" << nu->getWaitingTime() << "] | ";
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
void Interpolation<T>::displayFunction()
{
  cout << " - Core type : [ " << m_function->getCoreType() << " ]" << endl;
  cout << " - Cross section type : [ " << m_function->getCrossSection() << " ]" << endl;
}

template <typename T>
void Interpolation<T>::displayResults()
{
    double co_err_inf, tu_err_inf, ai_err_inf;
    double co_err_mse, tu_err_mse, ai_err_mse;
    co_err_inf = Utils::maxAbsValue(m_approxErrors[Cocagne]);
    tu_err_inf = Utils::maxAbsValue(m_approxErrors[Tucker]);
    ai_err_inf = Utils::maxAbsValue(m_approxErrors[AI]);
    co_err_mse = Utils::computeMseError(m_approxResults[Apollo],m_approxResults[Cocagne]);
    tu_err_mse = Utils::computeMseError(m_approxResults[Apollo],m_approxResults[Tucker]);
    ai_err_mse = Utils::computeMseError(m_approxResults[Tucker],m_approxResults[AI]);
    cout << " - Interpolation error using Cocagne (pcm)";
    cout << " [e_inf = " << co_err_inf << ", e_mse = " << co_err_mse << "]" << endl;
    cout << " - Interpolation error using Tucker (pcm)";
    cout << " [e_inf = " << tu_err_inf << ", e_mse = " << tu_err_mse << "]" << endl;
    cout << " - Interpolation error using AI (pcm)";
    cout << " [e_inf = " << ai_err_inf << ", e_mse = " << ai_err_mse << "]" << endl;
    Utils::separateur();
    cout << " - Number of calculation point = " << m_nbEvals << endl;
    cout << " - Total Time = " << m_totalTime << endl;
    cout << " - AI Run Time = " << m_runTime << endl;
}

template <typename T>
void Interpolation<T>::displayAll()
{
    Utils::separateur();
    Utils::separateur();
    cout << endl;
    displayFunction();
    displayRealDomain();
    Utils::separateur();
    displayInterpolationPoints();
    Utils::separateur();
    displayResults();
    cout << endl;
}
/******************************************************************************/

/*************************** Results ******************************************/
template <typename T>
void Interpolation<T>::computeAIApproximationResults()
{
    vector<double> ai_res, ai_err;
    double val, maxValue = Utils::maxAbsValue(m_approxResults[Tucker]);
    for (int j=0; j<m_nbTestPoints; j++)
    {
        val = interpolation(m_testPoints[j]);
        ai_res.push_back(val);
        ai_err.push_back(pow(10,5)*(val-m_approxResults[Tucker][j])/maxValue);
    }
    m_approxResults.insert(pair<method,vector<double>>(AI,ai_res));
    m_approxErrors.insert(pair<method,vector<double>>(AI,ai_err));
}

template <typename T>
void Interpolation<T>::readDataAndResults()
{
    readReferencePointsFromFile();
    readApproximationResultsFromFile();
}

template <typename T>
void Interpolation<T>::saveResults()
{
    computeAIApproximationResults();
    saveApproximationResultsInFile();
}

template <typename T>
void Interpolation<T>::readReferencePointsFromFile()
{
    string csName = m_function->getCrossSection();
    string s = Utils::replace(csName,"*","_");
    ifstream file(m_function->realDataDirPath() + "/FinalResults/" + s, ios::in);
    if(file)
    {
        string line;
        vector<double> data;
        vector<double> max(m_d,-numeric_limits<double>::max());
        vector<double> min(m_d,numeric_limits<double>::max());
        MultiVariatePoint<double> p(5,0);
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
    string csName = m_function->getCrossSection();
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
        m_approxResults.insert(pair<method,vector<double>>(Apollo,ap_res));
        m_approxResults.insert(pair<method,vector<double>>(Cocagne,co_res));
        m_approxResults.insert(pair<method,vector<double>>(Tucker,tu_res));
        m_approxErrors.insert(pair<method,vector<double>>(Cocagne,co_err));
        m_approxErrors.insert(pair<method,vector<double>>(Tucker,tu_err));
        file.close();
    }
    else cerr << "Error while opening " << csName << " file!" << endl;
}

template <typename T>
void Interpolation<T>::saveApproximationResultsInFile()
{
    string s = Utils::replace(m_function->getCrossSection(),"*","_");
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
            file << setprecision(m_precision) << m_approxResults[Apollo][i] << " ";
            file << setprecision(m_precision) << m_approxResults[Cocagne][i] << " ";
            file << setprecision(m_precision) << m_approxResults[Tucker][i] << " ";
            file << setprecision(m_precision) << m_approxResults[AI][i] << " ";
            file << setprecision(m_precision) << m_approxErrors[Cocagne][i] << " ";
            file << setprecision(m_precision) << m_approxErrors[Tucker][i] << " ";
            file << setprecision(m_precision) << m_approxErrors[AI][i] << endl;
        }
        file.close();
    }
    else cerr << "Error while opening the file!" << endl;
}
/******************************************************************************/

#endif
