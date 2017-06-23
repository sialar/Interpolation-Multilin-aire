#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/MixedInterpolation.hpp"
#include "../include/MultiVariatePoint.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    int dimD = (argc > 1) ? stoi(argv[1]) : 3;
    int dimN = (argc > 2) ? stoi(argv[2]) : 1;
    int nbTestPoints = (argc > 3) ? stoi(argv[3]) : 10;
    int maxIteration = (argc > 4) ? stoi(argv[4]) : 1000;
    MultiVariatePoint<double> p(dimD,0,0);
    MultiVariatePoint<int> methods(3,0,0);
    int f = (argc > 5) ? stoi(argv[5]) : 1;

    for (int i=0; i<dimD; i++)
        p(i) = (argc > 6+i) ? stoi(argv[6+i]) : Utils::randomValue(-1,1);
    for (int i=0; i<dimD; i++)
        methods(i) = (argc > 11+i) ? stoi(argv[11+i]) : Utils::randomValue(0,2);
    Function interpFunc;
    if (f==1) interpFunc = Functions::autoPolynomialFunction;
    if (f==2) interpFunc = Functions::functionToPlot;
    if (f==3) interpFunc = Functions::sinOfNorm2;
    if (f==4) interpFunc = Functions::h;
    if (f==5) interpFunc = Functions::phi;

    MixedInterpolationPtr interp(new MixedInterpolation(dimD,dimN,maxIteration,methods,interpFunc));
    interp->setSaveError(false);
    interp->disableProgressDisplay();

    // Initialisation of test points
    interp->setRandomTestPoints(nbTestPoints);

    // Path creation
    double execTime = interp->buildPathWithAIAlgo(chrono::steady_clock::now(), 1e-9, false);

    vector<vector<double>> realValues, estimate;
    for (MultiVariatePoint<double> p : interp->testPoints())
    {
        realValues.push_back(interp->func(p));
        estimate.push_back(interp->interpolation(p,interp->path().size()));
    }

    vector<double> res, exactValue;
    res = interp->interpolation(p,interp->path().size());
    exactValue = interpFunc(p,dimN);
    double relativeError = Utils::relativeInterpolationError(realValues,estimate);
    double mseError = Utils::mseInterpolationError(realValues,estimate);

    // Output
    ofstream file(Utils::projectPath + "data/res.txt", ios::out | ios::trunc);
    for (size_t i=0; i<res.size(); i++)
        file << res[i] << " ";
    file << endl;
    for (size_t i=0; i<exactValue.size(); i++)
        file << exactValue[i] << " ";
    file << endl;
    file << relativeError << endl;
    file << mseError << endl;
    file << execTime << endl;

    return 0;
}
