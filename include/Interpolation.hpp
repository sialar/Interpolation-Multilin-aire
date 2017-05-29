#ifndef INTERPOLATION
#define INTERPOLATION

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <memory>
#include <cmath>
#include <chrono>

#include "MultiVariatePoint.hpp"
#include "BinaryTree.hpp"
#include "Utils.hpp"

using namespace std;

class Interpolation
{
    protected:
        int m_d;
        int m_maxIteration;

        vector<vector<double>> m_interpolationPoints;
        vector<MultiVariatePoint<double>> m_interpolationNodes;
        vector<MultiVariatePoint<double>> m_testPoints;

        map<int, double> m_errors;
        bool m_saveError;

    public:
        Interpolation(int d, int nIter);
        virtual ~Interpolation() {};

        const vector<vector<double>>& points() { return m_interpolationPoints; };
        const vector<MultiVariatePoint<double>>& interpolationPoints() { return m_interpolationNodes; };
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };
        void setSaveError(bool error) { m_saveError = error; };
        bool getSaveError() { return m_saveError; };
        void displayInterpolationPointsInEachDirection();
        void displayInterpolationMultiVariatePoints();

        void saveErrorsInFile();
};

typedef std::unique_ptr<Interpolation> InterpolationPtr;

#endif
