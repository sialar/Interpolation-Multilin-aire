#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/PiecewiseInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int chooseDimension(int argc, char* argv[])
{
    if (argc > 1) return stoi(argv[1]);
    int dim = -1;
    while (dim < 0)
    {
        cout << " - Choose the space dimension : ";
        cin >> dim;
    }
    return dim;
}

int chooseNbTestPoints(int argc, char* argv[])
{
  if (argc > 2) return stoi(argv[2]);
  int nbTestPoints = -1;
  while (nbTestPoints < 0)
  {
    cout << " - Choose the number ot test points : ";
    cin >> nbTestPoints;
  }
  return nbTestPoints;
}

int chooseMethod(int argc, char* argv[])
{
    int method = -1;
    if (argc > 3) method = stoi(argv[3]);
    while (method!=1 && method!=2)
    {
        cout << " - Choose the method of interpolation: " << endl;
        cout << "\t - 1: Using piecewise functions and middle points: " << endl;
        cout << "\t - 2: Using quadratic functions and middle points: " << endl << " - ";
        cin >> method;
    }
    return method;
}

int chooseMaxIteration(int argc, char* argv[])
{
  if (argc > 4) return stoi(argv[4]);
  int maxIteration = -1;
  while (maxIteration < 0)
  {
    cout << " - Choose the maximum number of iteration : ";
    cin >> maxIteration;
  }
  Utils::separateur();
  return maxIteration;
}


int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dim = chooseDimension(argc,argv);
    int nbTestPoints = chooseNbTestPoints(argc,argv);
    int method = chooseMethod(argc,argv);
    int maxIteration = chooseMaxIteration(argc,argv);

    PiecewiseInterpolationPtr interp(new PiecewiseInterpolation(dim,maxIteration,method));

    // Initialisation of test points
    vector<MultiVariatePoint<double>> testPoints;
    vector<double> realValues, estimate;
    testPoints.resize(nbTestPoints);
    for (int j=0; j<nbTestPoints; j++)
        testPoints[j] = Utils::createRandomMultiVariatePoint(dim);
    interp->setTestPoints(testPoints);

    // Path creation
    double threshold = 1e-6;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = "
         << threshold;
    interp->testPathBuilt(threshold, maxIteration<21,0);

    // Computing real values, and approximation of function g at test points
    Utils::separateur();
    for (MultiVariatePoint<double> p : testPoints)
    {
        realValues.push_back(Utils::g(p));
        estimate.push_back(interp->interpolation_ND(p,interp->path().size()));
    }
    // Evaluation
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;

    /*
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = "
         << threshold;
    interp->testPathBuilt(threshold, maxIteration<21,1);
    Utils::separateur();
    for (MultiVariatePoint<double> p : testPoints)
    {
        realValues.push_back(Utils::f(p));
        estimate.push_back(interp->interpolation_ND(p,interp->path().size()));
    }
    // Evaluation
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;
    */
    Utils::separateur();
    return 0;
}
