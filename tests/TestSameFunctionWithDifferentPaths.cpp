#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/PiecewiseInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dimD = Utils::chooseDimensionD(argc,argv,1);
    int dimN = Utils::chooseDimensionN(argc,argv,2);
    int nbTestPoints = Utils::chooseNbTestPoints(argc,argv,3);
    int method = Utils::chooseMethod(argc,argv,4);
    int maxIteration = Utils::chooseMaxIteration(argc,argv,5);

    PiecewiseInterpolationPtr interp(new PiecewiseInterpolation(dimD,dimN,maxIteration,method,Utils::g));

    // Initialisation of test points
    interp->setRandomTestPoints(nbTestPoints);

    double threshold = 1e-6;
    cout << " - Interpolation of function g using the path of g" << endl;
    interp->testPathBuilt(threshold, maxIteration<21);
    interp->savePathInFile("data/path.txt");
    vector<vector<double>> realValues, estimate;
    for (MultiVariatePoint<double> p : interp->testPoints())
    {
        realValues.push_back(interp->func(p));
        estimate.push_back(interp->interpolation(p,interp->path().size()));
    }
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;

    Utils::separateur();
    interp->clearAllTrees();
    interp->clearAll();
    realValues.clear();
    estimate.clear();

    cout << " - Interpolation of function g using the path of f" << endl;
    cout << " - Computing the interpolation points obtained with function f" << endl;

    interp->setFunc(Utils::f);
    interp->testPathBuilt(threshold, maxIteration<21);
    interp->clearAllAlpha();
    interp->setFunc(Utils::g);
    interp->computeAllAlphaNuInPredefinedPath();
    interp->savePathInFile("data/other_path.txt");
    for (MultiVariatePoint<double> p : interp->testPoints())
    {
        realValues.push_back(interp->func(p));
        estimate.push_back(interp->interpolation(p,interp->path().size()));
    }
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;

    Utils::separateur();
    return 0;
}
