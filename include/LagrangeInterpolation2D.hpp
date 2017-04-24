#ifndef LAGRANGEINTERPOLATION2D
#define LAGRANGEINTERPOLATION2D

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <list>
#include <cmath>
#include <limits>
#include <time.h>

#include "Utils.hpp"
#include "IndiceND.hpp"

using namespace std;


class LagrangeInterpolation2D
{
    private:
        double m_alphaInitVal;
        vector<vector<double>> m_alphaTab;

        vector<double> m_pointsX;
        vector<double> m_pointsY;

        vector<IndiceND> m_path;

        list<IndiceND> m_curentNeighbours;

    public:
        LagrangeInterpolation2D(int sizeX, int sizeY,int path);
        ~LagrangeInterpolation2D();

        void clear();

        /************************* Data points ********************************/
        const vector<double>& pointsX() { return m_pointsX; };
        void setPointsX(vector<double> points) { m_pointsX = points; };
        const vector<double>& pointsY() { return m_pointsY; };
        void setPointsY(vector<double> points) { m_pointsY = points; };

        /************************* Path ***************************************/
        const vector<IndiceND>& path() { return m_path; };
        void choosePath(int n, int m, int v /* 0, 1 ou 2 */);
        int getLastIndiceInPath(int maxI, int maxJ);
        bool indiceInPath(IndiceND& index);
        double testPathBuilt(int nbIteration);
        void buildPathWithAIAlgo(int k);
        void savePathInFile();
        void showPath();

        /************************* Alpha **************************************/
        void setAlphaInitVal(double val) { m_alphaInitVal = val; };
        void initAlphaTab(double initVal);
        double computeOneAlphaNu(IndiceND& nu);
        double computeLastAlphaNu(IndiceND& nu);
        void computeAllAlphaNu();
        void showAlphaTab();

        /************************ Neighbours **********************************/
        void updateCurentNeighbours(IndiceND nu);
        bool isCorrectNeighbourToCurentPath(IndiceND nu);
        void showCurentNeighbours();

        /*********************** Interpolation ********************************/
        double lagrangeBasisFunction_1D(int j, int k, double y, int axis);
        double lagrangeInterpolation_2D_simple(double x, double y);
        double lagrangeInterpolation_2D_iterative(double x, double y);
        double computeExecTimeOfOneApprox();

};

#endif
