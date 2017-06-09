#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/MixedInterpolation.hpp"
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
    int dimD = Utils::chooseDimensionD(argc,argv,1);
    int dimN = Utils::chooseDimensionN(argc,argv,2);
    int nbTestPoints = Utils::chooseNbTestPoints(argc,argv,3);
    MultiVariatePoint<int> methods = chooseMethods(dimD);
    int maxIteration = Utils::chooseMaxIteration(argc,argv,4);
    MixedInterpolationPtr interp(new MixedInterpolation(dimD,dimN,maxIteration,methods,Utils::g));

    // Initialisation of test points
    interp->setRandomTestPoints(nbTestPoints);

    // Path creation
    double threshold = 1e-6;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = "
         << threshold;
    interp->testPathBuilt(threshold, maxIteration<21);

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
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;
    interp->savePathInFile("data/path.txt");
    Utils::separateur();

    return 0;
}
