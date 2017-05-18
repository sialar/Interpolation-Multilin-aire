#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/Interpolation.hpp"
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
  Utils::separateur();
  return nbTestPoints;
}

int chooseMethod(int argc, char* argv[])
{
    if (argc > 3) return stoi(argv[3]);
    int method = -1;
    while (method!=0 && method!=1)
    {
        cout << " - Choose the method of interpolation: " << endl;
        cout << "\t - 0: Using lagrange polynomial functions and Leja points " << endl;
        cout << "\t - 1: Using fine-tune functions and middle points: " << endl;
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

    Interpolation* interp = new Interpolation(dim,maxIteration,method);

    // Initialisation of test points
    vector<MultiVariatePoint<double>> testPoints;
    vector<double> realValues, estimate;
    testPoints.resize(nbTestPoints);
    for (int j=0; j<nbTestPoints; j++)
        testPoints[j] = Utils::createRandomMultiVariatePoint(dim);
    interp->setTestPoints(testPoints);
    cout << " - Sequence of " << nbTestPoints << " random test points : " ;
    Utils::displayPoints(testPoints);
    Utils::separateur();

    // Path creation
    double threshold = 1e-4;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = "
         << threshold << endl << endl;
    interp->testPathBuilt(threshold, false);
    Utils::separateur();

    // Computing real values at test points
    cout << " - Real values of function g evaluated at test points :" << endl;
    for (MultiVariatePoint<double> p : testPoints)
        realValues.push_back(Utils::gNd(p));
    Utils::displayPoints(realValues);

    // Approximating g at test points
    cout << endl << " - Approximation of function g at test points : " << endl;
    for (MultiVariatePoint<double> p : testPoints)
        estimate.push_back(interp->interpolation_ND(p));
    Utils::displayPoints(estimate);

    interp->storeInterpolationFunctions();
    interp->savePathInFile();
    interp->displayPath();

    Utils::separateur();
    interp->displayInterpolationMultiVariatePoints();
    cout << endl;
    interp->displayInterpolationPointsInEachDirection();
    Utils::separateur();

    // Evaluation
    cout << endl << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;
    Utils::separateur();

    return 0;
}
