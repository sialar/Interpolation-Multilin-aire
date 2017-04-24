#ifndef LAGRANGEINTERPOLATIONND
#define LAGRANGEINTERPOLATIONND

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <cmath>
#include "IndiceND.hpp"

using namespace std;

class LagrangeInterpolationND
{
    private:

        int m_d;

        double m_alphaInitVal;
        map<IndiceND, double> m_alphaTab;

        vector<vector<double>> m_points;

        set<IndiceND> m_path;

        list<IndiceND> m_curentNeighbours;

    public:
        LagrangeInterpolationND(vector<int> sizes);
        ~LagrangeInterpolationND();

        void clear();

        /************************* Data points ********************************/
        const vector<vector<double>>& points() { return m_points; };
        void setDirPoints(int i, vector<double> pointsI) { m_points[i] = pointsI; };
        void showPoints();

        /************************* Path ***************************************/
        const set<IndiceND>& path() { return m_path; };
        set<IndiceND>::iterator getLastIndiceInPath(vector<int> max);
        bool indiceInPath(IndiceND& index);
        void showPath();


        /************************* Alpha **************************************/
        //void setAlphaInitVal(double val) { m_alphaInitVal = val; };
        //void initAlphaTab(double initVal);
        //double computeLastAlphaNu(IndiceND& nu);

        /************************ Neighbours **********************************/
        //void updateCurentNeighbours(IndiceND nu);
        //bool isCorrectNeighbourToCurentPath(IndiceND nu);
        //void showCurentNeighbours();

        /*********************** Interpolation ********************************/
        //double lagrangeBasisFunction_1D(int j, int k, double y, int axis);
        //double lagrangeInterpolation_2D_iterative(double x, double y);
        //double computeExecTimeOfOneApprox();
};

#endif
