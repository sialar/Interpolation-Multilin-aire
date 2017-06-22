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
    int f = Utils::chooseFunction(argc,argv,6);
    bool save = Utils::saveResults(argc,argv,7);

    Function interpFunc;
    if (f==1) interpFunc = Functions::autoPolynomialFunction;
    if (f==2) interpFunc = Functions::functionToPlot;
    if (f==3) interpFunc = Functions::sinOfNorm2;
    if (f==4) interpFunc = Functions::h;
    if (f==5) interpFunc = Functions::phi;

    PiecewiseInterpolationPtr interp(new PiecewiseInterpolation(dimD,dimN,maxIteration,method,interpFunc));


    // Initialisation of test points
    interp->setRandomTestPoints(nbTestPoints);

    // Path creation
    double threshold = 1e-9;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = "
         << threshold;
    double execTime = interp->testPathBuilt(threshold, true);

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
    double relativeError = Utils::relativeInterpolationError(realValues,estimate);
    double mseError = Utils::mseInterpolationError(realValues,estimate);
    cout << " - Relative Interpolation error (pcm) = " << relativeError << endl;
    cout << " - MSE Interpolation error (pcm) = " << mseError << endl;
    cout << " - Number of evaluation = " << interp->nbEvals() << endl;

    if (save)
    {
        interp->saveInterpolationBasisFunctions();
        interp->saveInterpolationProgression();
        interp->savePathInFile(Utils::projectPath + "data/path.txt");
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


    ofstream file(Utils::projectPath + "data/res.txt", ios::out | ios::trunc);
    file << relativeError << endl;
    file << mseError << endl;
    file << execTime << endl;

    Utils::separateur();
    return 0;
}
