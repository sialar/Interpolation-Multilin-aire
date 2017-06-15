#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/LagrangeInterpolation.hpp"
#include "../include/PiecewiseInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    int dimD = (argc > 1) ? stoi(argv[1]) : 3;
    int dimN = (argc > 2) ? stoi(argv[2]) : 1;
    int nbTestPoints = (argc > 3) ? stoi(argv[3]) : 10;
    int method = (argc > 4) ? stoi(argv[4]) : 0;
    int maxIteration = (argc > 5) ? stoi(argv[5]) : 1000;
    MultiVariatePoint<double> p(dimD,0,0);
    p(0) = (argc > 6) ? stoi(argv[6]) : Utils::randomValue(-1,1);
    p(1) = (argc > 7) ? stoi(argv[7]) : Utils::randomValue(-1,1);
    p(2) = (argc > 8) ? stoi(argv[8]) : Utils::randomValue(-1,1);
    int f = (argc > 9) ? stoi(argv[9]) : 1;
    Function interpFunc;
    if (f==1) interpFunc = Functions::autoPolynomialFunction;
    if (f==2) interpFunc = Functions::functionToPlot;
    if (f==3) interpFunc = Functions::sinOfNorm2;

    LagrangeInterpolationPtr interp0(new LagrangeInterpolation(dimD,dimN,maxIteration,interpFunc));
    PiecewiseInterpolationPtr interp12(new PiecewiseInterpolation(dimD,dimN,maxIteration,method,interpFunc));

    interp0->disableProgressDisplay();
    interp12->disableProgressDisplay();

    // Initialisation of test points
    if (method) interp12->setRandomTestPoints(nbTestPoints);
    else interp0->setRandomTestPoints(nbTestPoints);

    double execTime;
    double threshold = 1e-10;
    if (method) execTime = interp12->testPathBuilt(threshold, maxIteration<21);
    else execTime = interp0->testPathBuilt(threshold, maxIteration<21);

    vector<vector<double>> realValues, estimate;
    for (int i=0; i<nbTestPoints; i++)
    {
        if (method) realValues.push_back(interp12->func(interp12->testPoints()[i]));
        else realValues.push_back(interp0->func(interp0->testPoints()[i]));
        if (method) estimate.push_back(interp12->interpolation(interp12->testPoints()[i],interp12->path().size()));
        else estimate.push_back(interp0->interpolation(interp0->testPoints()[i],interp0->path().size()));
    }

    vector<double> res, exactValue;
    if (method) res = interp12->interpolation(p,interp12->path().size());
    else res = interp0->interpolation(p,interp0->path().size());
    exactValue = interpFunc(p,dimN);
    double relativeError = Utils::relativeInterpolationError(realValues,estimate);
    double mseError = Utils::mseInterpolationError(realValues,estimate);

    if (method)
    {
        interp12->saveInterpolationBasisFunctions();
        interp12->saveInterpolationProgression();
        interp12->savePathInFile(Utils::projectPath + "data/path.txt");
    }
    else
    {
        interp0->saveInterpolationBasisFunctions();
        interp0->saveInterpolationProgression();
        interp0->savePathInFile(Utils::projectPath + "data/path.txt");
    }

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
