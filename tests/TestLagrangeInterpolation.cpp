#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/LagrangeInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dimD = Utils::chooseDimensionD(argc,argv,1);
    int dimN = Utils::chooseDimensionN(argc,argv,2);
    int nbTestPoints = Utils::chooseNbTestPoints(argc,argv,3);
    int maxIteration = Utils::chooseMaxIteration(argc,argv,4);
    bool save = Utils::saveResults(argc,argv,5);
    bool error = Utils::saveError(argc,argv,6);

    Function f = Functions::f;
    if (save) f = Functions::functionToPlot;
    LagrangeInterpolationPtr interp(new LagrangeInterpolation(dimD,dimN,maxIteration,f));
    interp->setSaveError(error);

    // Initialisation of test points
    interp->setRandomTestPoints(nbTestPoints);

    // Path creation
    double threshold = 1e-10;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = "
         << threshold;
    interp->testPathBuilt(threshold, maxIteration<11);

    // Computing real values, and approximation of function g at test points
    Utils::separateur();
    vector<vector<double>> realValues, estimate;
    for (MultiVariatePoint<double> p : interp->testPoints())
    {
        realValues.push_back(interp->func(p));
        estimate.push_back(interp->interpolation(p,interp->path().size()));
    }
    if (nbTestPoints < 11)
    {
        cout << " - Sequence of " << nbTestPoints << " random test points : " << endl;
        Utils::displayPoints(interp->testPoints());
        cout << " - Real values of function g evaluated at test points :" << endl;
        Utils::displayPoints(realValues);
        cout << " - Approximation of function g at test points : " << endl;
        Utils::displayPoints(estimate);
    }
    // Evaluation
    cout << " - Relative Interpolation error = " << Utils::relativeInterpolationError(realValues,estimate) << endl;
    cout << " - MSE Interpolation error = " << Utils::mseInterpolationError(realValues,estimate) << endl;

    if (save)
    {
        interp->saveInterpolationBasisFunctions();
        interp->saveInterpolationProgression();
        interp->savePathInFile("data/path.txt");
    }

    Utils::separateur();
    if (Utils::displayResults())
    {
      Utils::separateur();
      interp->displayPath();
      Utils::separateur();
      interp->displayInterpolationMultiVariatePoints();
      cout << endl;
      interp->displayInterpolationPointsInEachDirection();
    }

    Utils::separateur();
    return 0;
}
