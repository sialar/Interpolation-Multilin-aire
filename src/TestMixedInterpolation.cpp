#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/MixedInterpolation.hpp"
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

vector<int> chooseMethods(int dim)
{
    vector<int> methods(dim,-1);
    for (int i=0; i<dim; i++)
    {
        while (methods[i]!=0 && methods[i]!=1 && methods[i]!=2)
        {
            cout << " - Choose the method of interpolation in direction [" << i << "]: " << endl;
            cout << "\t - 0: Using lagrange polynomial functions and leja points: " << endl;
            cout << "\t - 1: Using piecewise functions and middle points: " << endl;
            cout << "\t - 2: Using quadratic functions and middle points: " << endl << " - ";
            cin >> methods[i];
        }
    }
    return methods;
}

int chooseMaxIteration(int argc, char* argv[])
{
  if (argc > 3) return stoi(argv[3]);
  int maxIteration = -1;
  while (maxIteration < 0)
  {
    cout << " - Choose the maximum number of iteration : ";
    cin >> maxIteration;
  }
  Utils::separateur();
  return maxIteration;
}

bool withBackup(int argc, char* argv[])
{
  if (argc > 4) return stoi(argv[4]);
  char store = 'x';
  while (store!='y' && store!='n')
  {
      cout << " - Store path and interpolation progression? (y/n) ";
      cin >> store;
  }
  return (store=='y');
}

bool saveError(int argc, char* argv[])
{
  if (argc > 5) return stoi(argv[5]);
  char e = 'x';
  while (e!='y' && e!='n')
  {
      cout << " - Store interpolation error at the end of the algorithm? (y/n) ";
      cin >> e;
  }
  return (e=='y');
}

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dim = chooseDimension(argc,argv);
    int nbTestPoints = chooseNbTestPoints(argc,argv);
    vector<int> methods = chooseMethods(dim);
    int maxIteration = chooseMaxIteration(argc,argv);
    bool store = withBackup(argc,argv);
    bool error = saveError(argc,argv);

    MixedInterpolationPtr interp(new MixedInterpolation(dim,maxIteration,methods));
    interp->setSaveError(error);

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
    interp->testPathBuilt(threshold, maxIteration<11);

    // Computing real values, and approximation of function g at test points
    Utils::separateur();
    for (MultiVariatePoint<double> p : testPoints)
    {
        realValues.push_back(Utils::gNd(p));
        estimate.push_back(interp->interpolation_ND(p,interp->path().size()));
    }
    if (nbTestPoints < 11)
    {
        cout << " - Sequence of " << nbTestPoints << " random test points : " << endl;
        Utils::displayPoints(testPoints);
        cout << " - Real values of function g evaluated at test points :" << endl;
        Utils::displayPoints(realValues);
        cout << " - Approximation of function g at test points : " << endl;
        Utils::displayPoints(estimate);
    }

    // Evaluation
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;

    if (store)
    {
        interp->storeInterpolationBasisFunctions();
        interp->storeInterpolationProgression();
        Utils::separateur();
        char plot = 'x';
        while (plot!='y' && plot!='n')
        {
          cout << " - Plot path: (y/n) " ;
          cin >> plot;
        }
        interp->savePathInFile(plot=='y');
    }

    Utils::separateur();
    char display = 'x';
    while (display!='y' && display!='n')
    {
      cout << " - Display path and interpolation points: (y/n) " ;
      cin >> display;
    }
    if (display=='y' /*maxIteration<pow(10,dim)+1*/)
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
