#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/LagrangeInterpolation.hpp"
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

bool withBackup()
{
  char store = 0;
  while (store!='y' && store!='n')
  {
      cout << " - Store path and interpolation progression ? (y/n) ";
      cin >> store;
  }
  return (store=='y');
}

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dim = chooseDimension(argc,argv);
    int nbTestPoints = chooseNbTestPoints(argc,argv);
    int maxIteration = chooseMaxIteration(argc,argv);

    LagrangeInterpolationPtr interp(new LagrangeInterpolation(dim,maxIteration));

    // Initialisation of test points
    vector<MultiVariatePoint<double>> testPoints;
    vector<double> realValues, estimate;
    testPoints.resize(nbTestPoints);
    for (int j=0; j<nbTestPoints; j++)
        testPoints[j] = Utils::createRandomMultiVariatePoint(dim);
    interp->setTestPoints(testPoints);

    // Path creation
    double threshold = 1e-4;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = "
         << threshold;
    interp->testPathBuilt(threshold, maxIteration<101);

    if (maxIteration<101)
    {
        interp->displayPath();
        Utils::separateur();
        interp->displayInterpolationMultiVariatePoints();
        cout << endl;
        interp->displayInterpolationPointsInEachDirection();
    }

    Utils::separateur();
    cout << " - Sequence of " << nbTestPoints << " random test points : " << endl;
    if (nbTestPoints < 11)
        Utils::displayPoints(testPoints);

    // Computing real values at test points
    cout << " - Real values of function g evaluated at test points :" << endl;
    for (MultiVariatePoint<double> p : testPoints)
        realValues.push_back(Utils::gNd(p));
    if (nbTestPoints < 11)
        Utils::displayPoints(realValues);

    // Approximating g at test points
    cout << " - Approximation of function g at test points : " << endl;
    for (MultiVariatePoint<double> p : testPoints)
        estimate.push_back(interp->interpolation_ND(p,interp->path().size()));
    if (nbTestPoints < 11) Utils::displayPoints(estimate);
    // Evaluation
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;
    Utils::separateur();

    if (withBackup())
    {
        interp->storeInterpolationBasisFunctions();
        interp->storeInterpolationProgression();
        interp->savePathInFile();
    }
    Utils::separateur();

    //cout << BinaryTree::getIndice(0.499999) << " " << BinaryTree::getValue(1835008) << endl;
    //cout << BinaryTree::getValue(3670016) << endl;

    return 0;
}
