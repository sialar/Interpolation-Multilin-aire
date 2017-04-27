#ifndef INTERPOLATIONND
#define INTERPOLATIONND

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <sys/time.h> // pour gettimeofday
#include <ctime> // Pour clock

#include "MultiVariatePoint.hpp"
#include "Utils.hpp"

#ifdef _OPENMP
   #include <pthread.h>
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
   #define omp_get_num_threads() 0
#endif

using namespace std;

class Interpolation
{
    private:

        int m_d;
        map<MultiVariatePoint<int>, double> m_alphaMap;
        vector<vector<double>> m_points;
        vector<MultiVariatePoint<int>> m_path;
        list<MultiVariatePoint<int>> m_curentNeighbours;

    public:

        Interpolation(vector<int> sizes);
        ~Interpolation();

        void clear();

        /************************* Data points ********************************/
        const vector<vector<double>>& points() { return m_points; };
        MultiVariatePoint<double> getPoint(MultiVariatePoint<int> nu);
        void setDirPoints(int i, vector<double> pointsI) { m_points[i] = pointsI; };
        void smartDiscretization();
        void displayPoints();

        /************************* Path ***************************************/
        const vector<MultiVariatePoint<int>>& path() { return m_path; };
        int getLastIndiceInPath(MultiVariatePoint<int> max);
        void buildPathWithAIAlgo(int k, bool parallel);
        double testPathBuilt(int nbIteration, bool parallel);
        bool indiceInPath(MultiVariatePoint<int>& index);
        void savePathInFile();
        void displayPath();

        /************************* Alpha **************************************/
        double computeLastAlphaNu(MultiVariatePoint<int>& nu);
        double computeLastAlphaNuPar(MultiVariatePoint<int>& nu);
        double testAlphaNuComputation(MultiVariatePoint<int>& nu);
        void displayAlphaTab();

        /************************ Neighbours **********************************/
        void updateCurentNeighbours(MultiVariatePoint<int>& nu);
        bool isCorrectNeighbourToCurentPath(MultiVariatePoint<int>& nu);
        void displayCurentNeighbours();

        /*********************** Interpolation ********************************/
        double lagrangeBasisFunction_1D(int j, int k, double t, int axis);
        double lagrangeInterpolation_ND_iterative(MultiVariatePoint<double> x);
        double computeExecTimeOfOneApprox();
        void displayInterpolationPoints();

};

#endif
