#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/LagrangeInterpolation.hpp"
#include "../include/PiecewiseInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;


int chooseMethod(int argc, char* argv[], int argNum)
{
    int method = -1;
    if (argc > 4) method = stoi(argv[4]);
    while (method!=0 && method!=1 && method!=2)
    {
        cout << " - Choose the method of interpolation: " << endl;
        cout << "\t - 0: Using Lagrange polynomial functions and Leja points: " << endl;
        cout << "\t - 1: Using piecewise functions and middle points: " << endl;
        cout << "\t - 2: Using quadratic functions and middle points: " << endl << " - ";
        cin >> method;
    }
    return method;
}

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dimD = Utils::chooseDimensionD(argc,argv,1);
    int dimN = Utils::chooseDimensionN(argc,argv,2);
    int nbTestPoints = Utils::chooseNbTestPoints(argc,argv,3);
    int method = chooseMethod(argc,argv,4);
    int maxIteration = Utils::chooseMaxIteration(argc,argv,5);
    int f = Utils::chooseFunction(argc,argv,6);

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

    double threshold = 1e-10;
    cout << " - Interpolation of function g using the path of g" << endl;
    if (method) interp12->testPathBuilt(threshold, maxIteration<21);
    else interp0->testPathBuilt(threshold, maxIteration<21);
    if (method) interp12->savePathInFile(Utils::projectPath + "data/path.txt");
    else interp0->savePathInFile(Utils::projectPath + "data/path.txt");

    vector<vector<double>> realValues, estimate;
    for (int i=0; i<nbTestPoints; i++)
    {
        if (method) realValues.push_back(interp12->func(interp12->testPoints()[i]));
        else realValues.push_back(interp0->func(interp0->testPoints()[i]));
        if (method) estimate.push_back(interp12->interpolation(interp12->testPoints()[i],interp12->path().size()));
        else estimate.push_back(interp0->interpolation(interp0->testPoints()[i],interp0->path().size()));
    }
    cout << " - Relative Interpolation error (pcm) = " << Utils::relativeInterpolationError(realValues,estimate) << endl;
    cout << " - MSE Interpolation error (pcm) = " << Utils::mseInterpolationError(realValues,estimate) << endl;

    Utils::separateur();
    if (method) interp12->clearAllTrees();
    if (method) interp12->clearAll();
    else interp0->clearAll();
    realValues.clear();
    estimate.clear();

    cout << " - Interpolation of function g using the path of f /= g" << endl;
    cout << " - Computing the interpolation points obtained with function f" << endl;

    if (method)
    {
      interp12->setFunc(Functions::h);
      interp12->testPathBuilt(threshold, maxIteration<21);
      interp12->clearAllAlpha();
      interp12->setFunc(interpFunc);
      interp12->computeAllAlphaNuInPredefinedPath();
      interp12->savePathInFile(Utils::projectPath + "data/other_path.txt");
    }
    else
    {
        interp0->setFunc(Functions::h);
        interp0->testPathBuilt(threshold, maxIteration<21);
        interp0->clearAllAlpha();
        interp0->setFunc(interpFunc);
        interp0->computeAllAlphaNuInPredefinedPath();
        interp0->savePathInFile(Utils::projectPath + "data/other_path.txt");
    }
    for (int i=0; i<nbTestPoints; i++)
    {
        if (method) realValues.push_back(interp12->func(interp12->testPoints()[i]));
        else realValues.push_back(interp0->func(interp0->testPoints()[i]));
        if (method) estimate.push_back(interp12->interpolation(interp12->testPoints()[i],interp12->path().size()));
        else estimate.push_back(interp0->interpolation(interp0->testPoints()[i],interp0->path().size()));
    }
    cout << " - Relative Interpolation error (pcm) = " << Utils::relativeInterpolationError(realValues,estimate) << endl;
    cout << " - MSE Interpolation error (pcm) = " << Utils::mseInterpolationError(realValues,estimate) << endl;

    Utils::separateur();
    return 0;
}
