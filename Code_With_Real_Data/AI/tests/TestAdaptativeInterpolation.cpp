#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/MixedInterpolation.hpp"
#include "../include/LagrangeInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;

MultiVariatePoint<int> chooseMethods(int dim)
{
    MultiVariatePoint<int> methods(dim,0,-1);
    for (int i=0; i<dim; i++)
    {
        while (methods(i)!=0 && methods(i)!=1 && methods(i)!=2)
        {
            cout << " - Choose the method of interpolation in direction [" << i << "]: " << endl;
            cout << "\t - 0: Using lagrange polynomial functions and leja points: " << endl;
            cout << "\t - 1: Using piecewise functions and middle points: " << endl;
            cout << "\t - 2: Using quadratic functions and middle points: " << endl << " - ";
            cin >> methods(i);
        }
    }
    return methods;
}

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dimD = 5;
    Function interpFunc = Utils::f;
    int maxIteration = Utils::chooseMaxIteration(argc,argv,1);
    LagrangeInterpolationPtr interp(new LagrangeInterpolation(dimD,1,maxIteration,interpFunc));
    interp->readEDFTestPointsFromFile();
    interp->displayRealDomain();
    Utils::separateur();

    // Path creation
    double threshold = 1e-9;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = " << threshold;
    double execTime = interp->testPathBuilt(threshold, false);

    // Computing real values, and approximation of function g at test points
    Utils::separateur();
    vector<vector<double>> realValues, estimate;
    for (MultiVariatePoint<double> p : interp->testPoints())
    {
        realValues.push_back(interp->func(p));
        estimate.push_back(interp->interpolation(p,interp->path().size()));
    }
    // Evaluation
    double relativeError = Utils::relativeInterpolationError(realValues,estimate);
    double mseError = Utils::mseInterpolationError(realValues,estimate);
    cout << " - Relative Interpolation error (pcm) = " << relativeError << endl;
    cout << " - MSE Interpolation error (pcm) = " << mseError << endl;
    cout << " - Number of evaluation = " << interp->nbEvals() << endl;
    cout << " - Run Time = " << execTime << endl;

    return 0;
}
