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

    PiecewiseInterpolationPtr interp(new PiecewiseInterpolation(dim,maxIteration,method,Utils::g));

    // Initialisation of test points
    vector<MultiVariatePoint<double>> testPoints;
    vector<double> realValues, estimate;
    testPoints.resize(nbTestPoints);
    for (int j=0; j<nbTestPoints; j++)
        testPoints[j] = Utils::createRandomMultiVariatePoint(dim);
    interp->setTestPoints(testPoints);

    double threshold = 1e-6;
    cout << " - Interpolation of function g using the path of g" << endl;
    interp->testPathBuilt(threshold, maxIteration<21);
    interp->savePathInFile("data/path.txt");
    for (MultiVariatePoint<double> p : testPoints)
    {
        realValues.push_back(interp->func(p));
        estimate.push_back(interp->interpolation_ND(p,interp->path().size()));
    }
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;


    Utils::separateur();
    interp->clearAllTrees();
    interp->clearAll();
    realValues.clear();
    estimate.clear();


    cout << " - Interpolation of function f using the path of g" << endl;
    cout << " - Computing the interpolation points obtained with function f" << endl;
    interp->testPathBuilt(threshold, maxIteration<21);
    interp->clearAllAlpha();
    interp->setFunc(Utils::f);
    interp->computeAllAlphaNuInPredefinedPath();
    interp->savePathInFile("data/other_path.txt");
    for (MultiVariatePoint<double> p : testPoints)
    {
        realValues.push_back(interp->func(p));
        estimate.push_back(interp->interpolation_ND(p,interp->path().size()));
    }
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;

    Utils::separateur();
    return 0;
}
