#ifndef INTERPOLATIONND
#define INTERPOLATIONND

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <cmath>
#include <sys/time.h>
#include <ctime>
#include <chrono>

#include "MultiVariatePoint.hpp"
#include "Utils.hpp"

using namespace std;

class Interpolation
{
    private:

        int m_method;
        int m_maxIteration;
        int m_d;

        vector<double> m_lejaSequence;
        vector<double> m_middlePoints;

        map<MultiVariatePoint<int>*, double> m_alphaMap;
        vector<set<double>> m_points;
        vector<MultiVariatePoint<double>> m_testPoints;
        vector<MultiVariatePoint<int>*> m_path;
        list<MultiVariatePoint<int>*> m_curentNeighbours;

        vector<BinaryTree*> m_middles;


    public:

        Interpolation(int d, int nIter, int method);
        ~Interpolation();

        /************************* Data points ********************************/
        const vector<set<double>>& points() { return m_points; };
        MultiVariatePoint<double> getPoint(MultiVariatePoint<int> nu);
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };
        void addInterpolationPoint(MultiVariatePoint<int> nu);

        /************************* AI algo ************************************/
        const vector<MultiVariatePoint<int>*>& path() { return m_path; };
        int buildPathWithAIAlgo(auto start_time, double threshold, bool debug);
        bool indiceInPath(MultiVariatePoint<int> index);
        double computeLastAlphaNu(MultiVariatePoint<int>* nu);
        void updateCurentNeighbours(MultiVariatePoint<int>* nu);
        void updateNextPoints(MultiVariatePoint<int>* nu);
        bool isCorrectNeighbourToCurentPath(MultiVariatePoint<int>* nu);

        /*********************** Interpolation ********************************/
        double piecewiseFunction_1D(int k, double t, int axis);
        double lagrangeBasisFunction_1D(int j, int k, double t, int axis);
        double interpolation_ND(MultiVariatePoint<double> x);
        void setMethod(int method) { m_method = method;};

        /********************** Test functions ********************************/
        void testPathBuilt(double threshold, bool debug);
        double tryWithCurentPath();

        /********************** Display functions *****************************/
        void displayPath();
        void displayAlphaTab();
        void displayCurentNeighbours();
        void displayInterpolationPoints();
        void savePathInFile();
        void storeInterpolationFunctions();
};

#endif
