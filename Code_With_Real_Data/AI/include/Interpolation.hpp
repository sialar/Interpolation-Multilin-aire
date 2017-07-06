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

        map<string,vector<double>> m_appoloResult;
        map<string,vector<double>> m_cocagneResult;
        map<string,vector<double>> m_tuckerResult;
        map<string,vector<double>> m_aiResult;
        map<string,vector<double>> m_cocagneError;
        map<string,vector<double>> m_tuckerError;
        map<string,vector<double>> m_aiError;

        vector<double> m_appoloReactivity;
        vector<double> m_cocagneReactivity;
        vector<double> m_tuckerReactivity;
        vector<double> m_aiReactivity;
        vector<double> m_cocagneReactivityError;
        vector<double> m_tuckerReactivityError;
        vector<double> m_aiReactivityError;

        bool m_displayProgress = true;
        FunctionsPtr m_function;

    public:

        Interpolation() {};
        virtual ~Interpolation() {};
        Interpolation(int d, string core, vector<string> cs, int nIter);

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

        void setFunc(string c);
        void setFunc(string c, vector<string> vr);
        FunctionsPtr getFunc() { return m_function; };
        vector<double> func(MultiVariatePoint<double> x);


        /************************* AI algo ************************************/
        void computeAIResults();
        void computeReactivity();
        void computeLastAlphaNu(MultiVariatePointPtr<T> nu);
        int buildPathWithAIAlgo( double threshold, bool debug);
        virtual MultiVariatePointPtr<T> maxElement(int iteration) = 0;
        virtual MultiVariatePointPtr<T> getFirstMultivariatePoint() = 0;
        virtual void updateCurentNeighbours(MultiVariatePointPtr<T> nu) = 0;
        virtual bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<T> nu) = 0;

        /************************* Interpolation ******************************/
        virtual double basisFunction_1D(T code, double t, int axis) = 0;
        vector<double> interpolation(MultiVariatePoint<double>& x, int end);
        vector<double> interpolation(MultiVariatePoint<double>& x);

        /************************* Other functions ****************************/
        void clearAll();
        void displayPath();
        void clearAllAlpha();
        void displayResults();
        void displayRealDomain();
        void saveReactivityInFile();
        void saveFinalResultsInFile();
        void readReactivityFromFile();
        void readTuckerDataFromFile();
        void displayCurentNeighbours();
        void displayCrossSectionNames();
        void displayInterpolationPoints();
};

template <typename T>
using InterpolationPtr = shared_ptr<Interpolation<T>>;

template <typename T>
Interpolation<T>::Interpolation(int d, string core, vector<string> cs, int nIter)
{
    m_d = d;
    m_n = cs.size();
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

template <typename T>
void Interpolation<T>::computeAIResults()
{
    for (int i=0; i<m_n; i++)
    {
        string csName = m_function->getCrossSections()[i];
        vector<double> ai_res, ai_err;
        double val, denom;
        denom = Utils::maxAbsValue(m_tuckerResult[csName]);
        for (int j=0; j<m_nbTestPoints; j++)
        {
            val = interpolation(m_testPoints[j])[i];
            ai_res.push_back(val);
            ai_err.push_back(pow(10,5) * (val-m_tuckerResult[csName][j]) / denom);
        }
        m_aiResult.insert(pair<string,vector<double>>(csName,ai_res));
        m_aiError.insert(pair<string,vector<double>>(csName,ai_err));
    }
    if (m_n == 12) computeReactivity();
}
template <typename T>
void Interpolation<T>::computeReactivity()
{
    double nu_f1, nu_f2, t1, t2, s011, s012, s021, s022;
    m_aiReactivity.resize(m_nbTestPoints);
    m_aiReactivityError.resize(m_nbTestPoints);
    map<string,vector<double>>::iterator it;
    for (int i=0; i<m_nbTestPoints; i++)
    {
        nu_f1 = m_aiResult["macro_nu*fission0"][i];
        nu_f2 = m_aiResult["macro_nu*fission1"][i];
        s011 = m_aiResult["macro_scattering000"][i];
        s012 = m_aiResult["macro_scattering001"][i];
        s021 = m_aiResult["macro_scattering010"][i];
        s022 = m_aiResult["macro_scattering011"][i];
        t1 = m_aiResult["macro_totale0"][i];
        t2 = m_aiResult["macro_totale1"][i];
        double num = nu_f1 * (t2 - s022) + nu_f2 * s012;
        double denom = (t1 - s011) * (t2 - s022) - s012 * s021;
        double k_eff = num / denom;
        m_aiReactivity[i] = 1 - 1/k_eff;
        m_aiReactivityError[i] = (m_aiReactivity[i]-m_tuckerReactivity[i]) * pow(10,5);
    }
    /*
    for (int i=0; i<m_nbTestPoints; i++)
    {
        nu_f1 = m_tuckerResult["macro_nu*fission0"][i];
        nu_f2 = m_tuckerResult["macro_nu*fission1"][i];
        s011 = m_tuckerResult["macro_scattering000"][i];
        s012 = m_tuckerResult["macro_scattering001"][i];
        s021 = m_tuckerResult["macro_scattering010"][i];
        s022 = m_tuckerResult["macro_scattering011"][i];
        t1 = m_tuckerResult["macro_totale0"][i];
        t2 = m_tuckerResult["macro_totale1"][i];
        double num = (nu_f1 * (t2 - s022)) + (nu_f2 * s012);
        double denom = ((t1 - s011) * (t2 - s022)) - (s012 * s021);
        double k_eff = num / denom;
        cout <<  (1-1/k_eff) << " " << m_tuckerReactivity[i] << endl;
        denom = Utils::maxAbsValue(m_appoloResult["macro_totale0"]);
        double e = pow(10,5) * (m_tuckerResult["macro_totale0"][i]-m_appoloResult["macro_totale0"][i]) / denom;
        cout << e << " " << m_tuckerError["macro_totale0"][i] << endl;
    }
    */
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
            cout << /*setprecision(numeric_limits<double>::digits10+1) <<*/ *it << " ";
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
    vector<double> co_err, tu_err, ai_err;
    for (int i=0; i<m_n; i++)
    {
        string csName = m_function->getCrossSections()[i];
        co_err.push_back(Utils::maxAbsValue(m_cocagneError[csName]));
        tu_err.push_back(Utils::maxAbsValue(m_tuckerError[csName]));
        ai_err.push_back(Utils::maxAbsValue(m_aiError[csName]));
    }

    cout << " - Interpolation error using Cocagne (pcm) = ";
    Utils::displayValues(co_err);
    cout << " - Interpolation error using Tucker (pcm) = ";
    Utils::displayValues(tu_err);
    cout << " - Interpolation error using AI (pcm) = ";
    Utils::displayValues(ai_err);
    double co_rea_err, tu_rea_err, ai_rea_err;
    if (m_n==12)
    {
        co_rea_err = Utils::maxAbsValue(m_cocagneReactivityError);
        tu_rea_err = Utils::maxAbsValue(m_tuckerReactivityError);
        ai_rea_err = Utils::maxAbsValue(m_aiReactivityError);
        cout << " - Reactivity error using Cocagne (pcm) = " << co_rea_err << endl;
        cout << " - Reactivity error using Tucker (pcm) = " << tu_rea_err << endl;
        cout << " - Reactivity error using AI (pcm) = " << ai_rea_err << endl;
    }
}

template <typename T>
void Interpolation<T>::readTuckerDataFromFile()
{
    for (int k=0; k<m_n; k++)
    {
        string csName = m_function->getCrossSections()[k];
        string s = Utils::replace(csName,"*","_");
        ifstream file(Utils::projectPath + "AI/data/" + m_function->getCoreType() + "/" + s, ios::in);
        if(file)
        {
            string line;
            vector<double> data, ap_res, co_res, tu_res, co_err, tu_err;
            vector<double> max(m_d,-numeric_limits<double>::max());
            vector<double> min(m_d,numeric_limits<double>::max());
            MultiVariatePoint<double> p(5,0,0);
            while (getline(file, line))
            {
                line = Utils::eraseExtraSpaces(line);
                stringstream ss(line);
                data.clear();
                string word;
                while (getline(ss, word, ' '))
                    data.push_back(stod(word));
                if (!k)
                {
                    for (int i=0; i<5; i++)
                        p(i) = data[i+1];
                    m_testPoints.push_back(p);
                    for (int i=0; i<m_d; i++)
                    {
                        if (p(i) > max[i])
                            max[i] = p(i);
                        if (p(i) < min[i])
                            min[i] = p(i);
                    }
                }
                ap_res.push_back(data[6]);
                co_res.push_back(data[7]);
                tu_res.push_back(data[8]);
                co_err.push_back(data[9]);
                tu_err.push_back(data[10]);
            }
            if (!k)
            {
                for (int i=0; i<m_d; i++)
                {
                    m_realDomain[i].push_back(min[i]);
                    m_realDomain[i].push_back(max[i]);
                }
                m_nbTestPoints = m_testPoints.size();
                for (int i=0; i<m_nbTestPoints; i++)
                    for (int j=0; j<m_d; j++)
                        m_testPoints[i](j) = Utils::adaptCoordsToFunctionDomain(min[j], max[j], m_testPoints[i](j));
            }
            m_appoloResult.insert(pair<string,vector<double>>(csName,ap_res));
            m_cocagneResult.insert(pair<string,vector<double>>(csName,co_res));
            m_tuckerResult.insert(pair<string,vector<double>>(csName,tu_res));
            m_cocagneError.insert(pair<string,vector<double>>(csName,co_err));
            m_tuckerError.insert(pair<string,vector<double>>(csName,tu_err));
            file.close();
        }
        else
            cerr << "Error while opening " << csName << " file!" << endl;
    }
    if (m_n==12) readReactivityFromFile();
}

template <typename T>
void Interpolation<T>::saveFinalResultsInFile()
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

              file << setprecision(numeric_limits<double>::digits10+1) << (i+1) << " ";
              for (int j=0; j<m_d; j++)
                  file << setprecision(numeric_limits<double>::digits10+1)  << x(j) << " ";
              file << setprecision(numeric_limits<double>::digits10+1)  << m_appoloResult[csName][i] << " " << m_cocagneResult[csName][i] << " " << m_tuckerResult[csName][i] << " ";
              file << setprecision(numeric_limits<double>::digits10+1)  << m_cocagneError[csName][i] << " " << m_tuckerError[csName][i] << " " << m_aiResult[csName][i] << " " << m_aiError[csName][i] << endl;//" " << Utils::maxAbsValue(m_appoloResult[csName]) << " " << Utils::maxAbsValue(m_tuckerResult[csName]) << " " << Utils::maxAbsValue(m_aiResult[csName]) << endl;
            }
            file.close();
        }
        else
            cerr << "Error while opening the file!" << endl;
    }
    if (m_n==12) saveReactivityInFile();
}

template <typename T>
void Interpolation<T>::readReactivityFromFile()
{
    ifstream file(Utils::projectPath + "AI/data/" + m_function->getCoreType() + "/ReactivityError", ios::in );
    if(file)
    {
        string line;
        vector<double> data;
        while (getline(file, line))
        {
            line = Utils::eraseExtraSpaces(line);
            stringstream ss(line);
            data.clear();
            string word;
            while (getline(ss, word, ' '))
                data.push_back(stod(word));
            m_appoloReactivity.push_back(data[6]);
            m_cocagneReactivity.push_back(data[7]);
            m_tuckerReactivity.push_back(data[8]);
            m_cocagneReactivityError.push_back(data[9]);
            m_tuckerReactivityError.push_back(data[10]);
        }
        file.close();
    }
    else
        cerr << "Error while opening the file!" << endl;
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

          file << setprecision(numeric_limits<double>::digits10+1) << (i+1) << " ";
          for (int j=0; j<m_d; j++)
              file << setprecision(numeric_limits<double>::digits10+1) << x(j) << " ";
          file << setprecision(numeric_limits<double>::digits10+1) << m_appoloReactivity[i] << " " << m_cocagneReactivity[i] << " " << m_tuckerReactivity[i] << " ";
          file << setprecision(numeric_limits<double>::digits10+1) << m_cocagneReactivityError[i] << " " << m_tuckerReactivityError[i] << " " << m_aiReactivity[i] << " ";
          file << setprecision(numeric_limits<double>::digits10+1) << m_aiReactivityError[i] << endl;
        }
        file.close();
    }
    else
        cerr << "Error while opening the file!" << endl;
}
/******************************************************************************/

#endif
