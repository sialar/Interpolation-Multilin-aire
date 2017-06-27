#ifndef INTERPOLATION
#define INTERPOLATION

#include <iostream>
#include <vector>
#include <list>
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

        bool m_displayProgress = true;

        FunctionsPtr m_function;

    public:
        Interpolation() {};
        virtual ~Interpolation() {};
        Interpolation(int d, int m_n, int nIter);

        /************************* Data points ********************************/
        const int nbEvals() { return m_nbEvals; };
        const double runTime() { return m_runTime; };
        const double totalTime() { return m_totalTime; };
        const int maxIteration() { return m_maxIteration; };
        void disableProgressDisplay() { m_displayProgress = false; };
        const vector<MultiVariatePointPtr<T>>& path() { return m_path; };
        virtual void addInterpolationPoint(MultiVariatePoint<double> p) = 0;
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        virtual MultiVariatePoint<double> getPoint(MultiVariatePointPtr<T> nu) = 0;
        const vector<MultiVariatePoint<double>>& interpolationPoints() { return m_interpolationNodes; };

        void setRandomTestPoints(int nbTestPoints);
        vector<MultiVariatePoint<double>> testPoints() { return m_testPoints; };
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };

        void setFunc(string c);
        void setFunc(string c, vector<string> vr);
        vector<double> func(MultiVariatePoint<double> x);


        /************************* AI algo ************************************/
        void computeLastAlphaNu(MultiVariatePointPtr<T> nu);
        int buildPathWithAIAlgo( double threshold, bool debug);
        virtual MultiVariatePointPtr<T> maxElement(int iteration) = 0;
        virtual MultiVariatePointPtr<T> getFirstMultivariatePoint() = 0;
        virtual void updateCurentNeighbours(MultiVariatePointPtr<T> nu) = 0;
        virtual bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<T> nu) = 0;

        /************************* Interpolation ******************************/
        virtual double basisFunction_1D(T code, double t, int axis) = 0;
        vector<double> interpolation(MultiVariatePoint<double>& x, int end);

        /************************* Other functions ****************************/
        void clearAll();
        void displayPath();
        void clearAllAlpha();
        void displayRealDomain();
        void displayCurentNeighbours();
        void readEDFTestPointsFromFile();
};

template <typename T>
using InterpolationPtr = shared_ptr<Interpolation<T>>;

template <typename T>
Interpolation<T>::Interpolation(int d, int n, int nIter)
{
    m_interpolationPoints.resize(d);
    m_function = make_shared<Functions>();
    m_realDomain.resize(d);
    m_maxIteration = nIter;
    m_d = d;
    m_n = n;
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
    m_function->setCoreType(c);
    m_function->setReactionTypes(vr);
}

template <typename T>
void Interpolation<T>::setFunc(string c)
{
    m_function->setCoreType(c);
    m_function->setAllReactionTypes();
}

template <typename T>
vector<double> Interpolation<T>::func(MultiVariatePoint<double> x)
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
int Interpolation<T>::buildPathWithAIAlgo(double threshold, bool debug)
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
/******************************************************************************/

/*************************** Display function *********************************/

template <typename T>
void Interpolation<T>::displayPath()
{
    // Format: (nu:points[nu]) --> () ...
    cout << " - Path =";
    int n = m_path.size();
    for (int i=0; i<n; i++)
    {
        if (i>0) cout << "\t";
        cout << " " << i << " :";
        cout << " [" << *m_path[i] << ":" << setprecision(numeric_limits<double>::digits10+1) << getPoint(m_path[i]);
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
        cout << "(" << (*nu) << ":" << setprecision(numeric_limits<double>::digits10+1) << getPoint(nu) << ":";
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
void Interpolation<T>::readEDFTestPointsFromFile()
{
    ifstream file(Utils::projectPath + "AI/data/reference_points.dat", ios::in);
    if(file)
    {
        string line;
        vector<double> max(m_d,-numeric_limits<double>::max());
        vector<double> min(m_d,numeric_limits<double>::max());
        MultiVariatePoint<double> p;
        while (getline(file, line))
        {
            p = Utils::getCoordsFromString(line);
            m_testPoints.push_back(p);
            for (int i=0; i<m_d; i++)
            {
                if (p(i) > max[i])
                    max[i] = p(i);
                if (p(i) < min[i])
                    min[i] = p(i);
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
    else
        cerr << "Error while opening the file!" << endl;
}
/******************************************************************************/

#endif
