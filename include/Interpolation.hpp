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
        int nbMethods = 3;

        vector<vector<double>> m_interpolationPoints;
        vector<MultiVariatePoint<double>> m_interpolationNodes;
        vector<MultiVariatePoint<double>> m_testPoints;

        vector<MultiVariatePointPtr<T>> m_path;
        list<MultiVariatePointPtr<T>> m_curentNeighbours;

        map<int, double> m_errors;
        bool m_displayProgress = true;
        bool m_saveError;

        Function m_function;

    public:
        Interpolation() {};
        virtual ~Interpolation() {};
        Interpolation(int d, int m_n, int nIter, Function f);

        /************************* Data points ********************************/
        const int maxIteration() { return m_maxIteration; };
        const vector<MultiVariatePoint<double>>& interpolationPoints() { return m_interpolationNodes; };
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        virtual MultiVariatePoint<double> getPoint(MultiVariatePointPtr<T> nu) = 0;
        virtual void addInterpolationPoint(MultiVariatePoint<double> p) = 0;

        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };
        vector<MultiVariatePoint<double>> testPoints() { return m_testPoints; };
        void setRandomTestPoints(int nbTestPoints);

        vector<double> func(MultiVariatePoint<double> x) { return m_function(x,m_n); };
        void setFunc(Function f) { m_function = f; };
        void setSaveError(bool error) { m_saveError = error; };
        map<int, double>& errors() { return m_errors; };
        void disableProgressDisplay() { m_displayProgress = false; };

        /************************* AI algo ************************************/
        vector<double> tryWithCurentPath();
        const vector<MultiVariatePointPtr<T>>& path() { return m_path; };
        double testPathBuilt(double threshold, bool debug);
        int buildPathWithAIAlgo(auto start_time, double threshold, bool debug);
        void computeLastAlphaNu(MultiVariatePointPtr<T> nu);
        void computeAllAlphaNuInPredefinedPath();
        virtual MultiVariatePointPtr<T> getFirstMultivariatePoint() = 0;
        virtual MultiVariatePointPtr<T> maxElement(int iteration) = 0;
        virtual void updateCurentNeighbours(MultiVariatePointPtr<T> nu) = 0;
        virtual bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<T> nu) = 0;

        /************************* Interpolation ******************************/
        vector<double> interpolation(MultiVariatePoint<double>& x, int end);
        virtual double basisFunction_1D(T code, double t, int axis) = 0;

        /************************* Display function ***************************/
        void displayInterpolationPointsInEachDirection();
        void displayInterpolationMultiVariatePoints();
        void savePathInFile(string fileName);
        void saveInterpolationProgression();
        void displayCurentNeighbours();
        void saveTestPointsInFile();
        void saveErrorsInFile();
        void clearAllAlpha();
        void displayPath();
        void displayAll();
        void clearAll();
};

template <typename T>
using InterpolationPtr = shared_ptr<Interpolation<T>>;

template <typename T>
Interpolation<T>::Interpolation(int d, int n, int nIter, Function f)
{
    m_interpolationPoints.resize(d);
    m_maxIteration = nIter;
    Functions::setCoefs(7,d,n);
    m_saveError = false;
    m_function = f;
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

/*************************** Data Points **************************************/
template <typename T>
void Interpolation<T>::setRandomTestPoints(int nbTestPoints)
{
  vector<MultiVariatePoint<double>> testPoints;
  testPoints.resize(nbTestPoints);
  for (int j=0; j<nbTestPoints; j++)
      testPoints[j] = Utils::createRandomMultiVariatePoint(m_d);
  setTestPoints(testPoints);
  //saveTestPointsInFile();
}
/******************************************************************************/

/***************************** AI algo ****************************************/
template <typename T>
int Interpolation<T>::buildPathWithAIAlgo(auto start_time, double threshold, bool debug)
{
    debug = false;
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

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI
        int fact = (m_saveError) ? 100 : 10;
        int step = floor(m_maxIteration/fact);
        if (!step) step = 1;
        if (iteration%step==0)
        {
            auto end_time = chrono::steady_clock::now();
            std::chrono::duration<double> run_time = end_time - start_time;
            vector<double> errors = tryWithCurentPath();
            if (m_saveError) m_errors.insert(pair<int, double>(iteration, errors[0]));
            if (m_displayProgress)
            {
                cout << endl << "\t- Interpolation error after " << iteration << " iterations: ";
                cout << "(Relative_e = " << errors[0]; //<< ", MSE_e = " << errors[1];
                cout << ") | Elapsed time : "  << run_time.count();
                if (errors[0] < threshold /*&& errors[1] < threshold*/)
                {
                    cout << endl << "   - AI Algo stoped after " << iteration << " iterations";
                    cout << " | Elapsed time : "  << run_time.count() << endl;
                    return iteration;
                }
            }
            saveErrorsInFile();
        }
        iteration++;
    }
    return iteration;
}
template <typename T>
void Interpolation<T>::computeLastAlphaNu(MultiVariatePointPtr<T> nu)
{
    double basisFuncProd = 1.0;
    vector<double> res = func(getPoint(nu));
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
void Interpolation<T>::computeAllAlphaNuInPredefinedPath()
{
    double basisFuncProd = 1.0;
    vector<double> res(m_n,0.0);
    for (int i=0; i<int(m_path.size()); i++)
    {
        res = func(getPoint(m_path[i]));
        for (int j=0; j<i; j++)
        {
            basisFuncProd = 1.1;
            for (int p=0; p<m_d; p++)
                basisFuncProd *= basisFunction_1D((*m_path[j])(p),getPoint(m_path[i])(p),p);
            for (int k=0; k<m_n; k++)
                res[k] -= m_path[j]->getAlpha()[k] * basisFuncProd;
        }
        m_path[i]->setAlpha(res);
    }
}
template <typename T>
double Interpolation<T>::testPathBuilt(double threshold, bool debug)
{
  auto start_time = chrono::steady_clock::now();
  int nbIterations = buildPathWithAIAlgo(start_time, threshold, debug);
  auto end_time = chrono::steady_clock::now();
  std::chrono::duration<double> run_time = end_time - start_time;
  cout << endl << "   - Time required to compute the path with " << nbIterations <<
  " iterations: " << run_time.count() << "s" << endl;
  return run_time.count();
}
template <typename T>
vector<double> Interpolation<T>::tryWithCurentPath()
{
  vector<vector<double>> realValues, estimate;
  vector<double> errors;
  for (MultiVariatePoint<double> p : m_testPoints)
  {
    realValues.push_back(func(p));
    estimate.push_back(interpolation(p,m_path.size()));
  }
  errors.push_back(Utils::relativeInterpolationError(realValues,estimate));
  errors.push_back(Utils::mseInterpolationError(realValues,estimate));
  return errors;
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
void Interpolation<T>::displayAll()
{
    displayPath();
    displayCurentNeighbours();
    displayInterpolationMultiVariatePoints();
    displayInterpolationPointsInEachDirection();
}

template <typename T>
void Interpolation<T>::saveInterpolationProgression()
{
  ofstream file(Utils::projectPath + "data/interpolation_progression.txt", ios::out | ios::trunc);
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
void Interpolation<T>::saveTestPointsInFile()
{
    ofstream file(Utils::projectPath + "data/test_points.txt", ios::out | ios::trunc);
    if(file)
    {
        for (MultiVariatePoint<double> x : m_testPoints)
        {
            for (int i=0; i<m_d-1; i++)
                file << x(i) << " ";
            file << x(m_d-1) << endl;
        }
        file.close();
    }
    else
        cerr << "Error while opening the file!" << endl;
}

template <typename T>
void Interpolation<T>::saveErrorsInFile()
{
  ofstream file(Utils::projectPath + "data/interpolation_error.txt", ios::out | ios::trunc);
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
void Interpolation<T>::savePathInFile(string fileName)
{
  ofstream file(fileName, ios::out | ios::trunc);
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
            for (int i=0; i<m_d; i++)
                file << (*nu)(i) << " " ;
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
/******************************************************************************/


#endif
