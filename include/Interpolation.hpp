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
#include "Utils.hpp"

using namespace std;

template <typename T>
class Interpolation
{
    protected:
        int m_d;
        int m_maxIteration;
        int nbMethods = 3;

        vector<vector<double>> m_interpolationPoints;
        vector<MultiVariatePoint<double>> m_interpolationNodes;
        vector<MultiVariatePoint<double>> m_testPoints;

        vector<MultiVariatePointPtr<T>> m_path;
        map<MultiVariatePointPtr<T>, double> m_alphaMap;
        list<MultiVariatePointPtr<T>> m_curentNeighbours;

        map<int, double> m_errors;
        bool m_saveError;

    public:
        Interpolation(int d, int nIter);
        virtual ~Interpolation() {};

        /************************* Data points ********************************/
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        const vector<MultiVariatePoint<double>>& interpolationPoints() { return m_interpolationNodes; };
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };
        void setSaveError(bool error) { m_saveError = error; };
        virtual MultiVariatePoint<double> getPoint(MultiVariatePointPtr<T> nu) = 0;
        virtual void addInterpolationPoint(MultiVariatePoint<double> p) = 0;

        /************************* AI algo ************************************/
        double tryWithCurentPath(int function);
        const vector<MultiVariatePointPtr<T>>& path() { return m_path; };
        void testPathBuilt(double threshold, bool debug, int function);
        int buildPathWithAIAlgo(auto start_time, double threshold, bool debug, int function);
        double computeLastAlphaNu(MultiVariatePointPtr<T> nu, int function);
        virtual MultiVariatePointPtr<T> getFirstMultivariatePoint() = 0;
        virtual MultiVariatePointPtr<T> maxElement(int iteration) = 0;
        virtual void updateCurentNeighbours(MultiVariatePointPtr<T> nu) = 0;
        virtual bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<T> nu) = 0;

        /************************* Interpolation ******************************/
        double interpolation_ND(MultiVariatePoint<double>& x, int end);
        virtual double basisFunction_1D(T code, double t, int axis) = 0;

        /************************* Display function ***************************/
        void displayInterpolationPointsInEachDirection();
        void displayInterpolationMultiVariatePoints();
        void storeInterpolationProgression();
        void displayCurentNeighbours();
        void saveErrorsInFile();
        void savePathInFile();
        void displayPath();
        void displayAll();
        void clearAll();
};

template <typename T>
using InterpolationPtr = shared_ptr<Interpolation<T>>;

template <typename T>
Interpolation<T>::Interpolation(int d, int nIter) : m_d(d), m_maxIteration(nIter)
{
    m_interpolationPoints.resize(m_d);
    m_saveError = false;
}

template <typename T>
void Interpolation<T>::clearAll()
{
    m_path.clear();
    m_alphaMap.clear();
    m_curentNeighbours.clear();
    m_interpolationNodes.clear();
    for (int i=0; i<m_d; i++)
        m_interpolationPoints[i].clear();
    m_interpolationPoints.clear();
}
/***************************** AI algo ****************************************/
template <typename T>
int Interpolation<T>::buildPathWithAIAlgo(auto start_time, double threshold, bool debug, int function)
{
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
                computeLastAlphaNu(nu, function);

        if (debug)
        {
            Utils::separateur();
            displayPath();
            displayCurentNeighbours();
        }

        argmax = maxElement(iteration);
        m_path.push_back(argmax);
        addInterpolationPoint(getPoint(argmax));
        updateCurentNeighbours(argmax);
        iteration++;

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI
        int fact = (m_saveError) ? m_maxIteration : 10;
        if (m_maxIteration>=10 && iteration%(m_maxIteration/fact)==0)
        {
            auto end_time = chrono::steady_clock::now();
            std::chrono::duration<double> run_time = end_time - start_time;
            double error = tryWithCurentPath(function);
            if (m_saveError) m_errors.insert(pair<int, double>(iteration, error));
            cout << endl << "\t- Interpolation error after " << iteration << " iterations: " << error;
            cout << " | Elapsed time : "  << run_time.count();
            if (error < threshold)
            {
                cout << "   - AI Algo stoped after " << iteration << " iterations";
                cout << " | Elapsed time : "  << run_time.count() << endl;
                return iteration;
            }
            saveErrorsInFile();
        }
    }
    return iteration;
}
template <typename T>
double Interpolation<T>::computeLastAlphaNu(MultiVariatePointPtr<T> nu, int function)
{
    double basisFuncProd;
    double res = 0;
    if (function) res = Utils::f(getPoint(nu));
    else res = Utils::g(getPoint(nu));
    for (MultiVariatePointPtr<T> l : m_path)
    {
        basisFuncProd = 1;
        for (int p=0; p<m_d; p++)
            basisFuncProd *= basisFunction_1D((*l)(p),getPoint(nu)(p),p);
        res -= l->getAlpha() * basisFuncProd;
    }
    m_alphaMap.insert(pair<MultiVariatePointPtr<T>,double>(nu,res));
    nu->setAlpha(res);
    return res;
}
template <typename T>
void Interpolation<T>::testPathBuilt(double threshold, bool debug, int function)
{
  auto start_time = chrono::steady_clock::now();
  int nbIterations = buildPathWithAIAlgo(start_time, threshold, debug, function);
  auto end_time = chrono::steady_clock::now();
  std::chrono::duration<double> run_time = end_time - start_time;
  cout << endl << "   - Time required to compute the path with " << nbIterations <<
  " iterations: " << run_time.count() << "s" << endl;
}
template <typename T>
double Interpolation<T>::tryWithCurentPath(int function)
{
  vector<double> realValues, estimate;
  for (MultiVariatePoint<double> p : m_testPoints)
  {
    if (function) realValues.push_back(Utils::f(p));
    else realValues.push_back(Utils::g(p));
    estimate.push_back(interpolation_ND(p,m_path.size()));
  }
  return Utils::interpolationError(realValues,estimate);
}
/******************************************************************************/

/*************************** Interpolation ************************************/
template <typename T>
double Interpolation<T>::interpolation_ND(MultiVariatePoint<double>& x, int end)
{
  double l_prod, sum = 0;
  for (int k=0; k<end; k++)
  {
      l_prod = 1;
      for (int i=0; i<m_d; i++)
          l_prod *= basisFunction_1D((*m_path[k])(i),x(i),i);
      sum += l_prod * m_path[k]->getAlpha();
  }
  return sum;
}
/******************************************************************************/

/*************************** Display function *********************************/
template <typename T>
void Interpolation<T>::displayInterpolationPointsInEachDirection()
{
    vector<double>::iterator it;
    for (int i=0; i<m_d; ++i)
    {
        cout << " - " << m_interpolationPoints[i].size() << " points in direction " << i << " : { ";
        for (it=m_interpolationPoints[i].begin(); it!=m_interpolationPoints[i].end(); it++)
            cout << setprecision(numeric_limits<double>::digits10+1) << *it << " ";
        cout << "}" << endl;
    }
}
template <typename T>
void Interpolation<T>::displayInterpolationMultiVariatePoints()
{
    cout << " - Interpolation nodes: { ";
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
        if (i>0) cout << "\t";
        cout << " " << i << " :";
        cout << " [" << *m_path[i] << ":" << setprecision(numeric_limits<double>::digits10+1) << getPoint(m_path[i]);
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
        cout << "(" << (*nu) << ":" << setprecision(numeric_limits<double>::digits10+1) << getPoint(nu) << ":";
        cout << nu->getAlpha() << ":" << nu << ") [" << nu->getWaitingTime() << "] | ";
    }
    cout << endl << endl;
}

template <typename T>
void Interpolation<T>::displayAll()
{
    displayPath();
    displayCurentNeighbours();
    displayInterpolationMultiVariatePoints();
    displayInterpolationPointsInEachDirection();
}

template <typename T>
void Interpolation<T>::storeInterpolationProgression()
{
  ofstream file("data/interpolation_progression.txt", ios::out | ios::trunc);
  if(file)
  {
      if (m_d==1)
      {
          vector<double> x;
          int nbPoints = int(m_testPoints.size());
          for (int i=0; i<nbPoints; i++)
              x.push_back(m_testPoints[i](0));
          sort(x.begin(), x.end());
          MultiVariatePoint<double> p;
          vector<double> tempPath;
          for (int j=0; j<nbPoints; j++)
          {
              for (int i=0; i<int(m_path.size()); i++)
              {
                  p = MultiVariatePoint<double>::toMonoVariatePoint(x[j]);
                  file << interpolation_ND(p, i+1) << " ";
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
void Interpolation<T>::saveErrorsInFile()
{
  ofstream file("data/interpolation_error.txt", ios::out | ios::trunc);
  if(file)
  {
      file << m_errors.size() << endl;
      map<int,double>::iterator it;
      for (it=m_errors.begin(); it!=m_errors.end(); it++)
          file << get<0>(*it) << " " << get<1>(*it) << endl;
      file.close();
  }
  else
      cerr << "Error while opening the file!" << endl;
}
template <typename T>
void Interpolation<T>::savePathInFile()
{
  ofstream file("data/path.txt", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==2 || m_d==3)
    {
        for (int i=0; i<m_d; i++)
            file << m_interpolationPoints[i].size() << " ";
        file << endl << m_path.size() << endl;
        MultiVariatePoint<double> x(m_d,0.0);
        for (MultiVariatePointPtr<T> nu : m_path)
        {
            x = getPoint(nu);
            for (int i=0; i<m_d; i++)
                file << x(i) << " ";
            for (int i=0; i<m_d; i++)
                file << (*nu)(i) << " " ;
            file << endl;
        }
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}
/******************************************************************************/


#endif
