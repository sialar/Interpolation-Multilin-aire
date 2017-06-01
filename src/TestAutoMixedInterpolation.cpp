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

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dim = chooseDimension(argc,argv);
    int nbTestPoints = chooseNbTestPoints(argc,argv);
    MultiVariatePoint<int> methods(dim,0);
    int maxIteration = chooseMaxIteration(argc,argv);

    MixedInterpolationPtr interp(new MixedInterpolation(dim,maxIteration,methods));

    // Initialisation of test points
    vector<MultiVariatePoint<double>> testPoints;
    vector<double> realValues, estimate;
    testPoints.resize(nbTestPoints);
    for (int j=0; j<nbTestPoints; j++)
        testPoints[j] = Utils::createRandomMultiVariatePoint(dim);
    interp->setTestPoints(testPoints);

    // Path creation
    double threshold = 1e-10;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will be tested with different methods in each direction." << endl;
    cout << " - Each time, the algorithm will stop when the interpolation error becomes lower than"
         << " a threshold = " << threshold << endl;

    MultiVariatePoint<int> optimalMethods = interp->tryAllCases(threshold,0);
    interp->setMethods(optimalMethods);

    Utils::separateur();
    interp->tryWithDifferentMethods(optimalMethods, threshold,0);
    interp->savePathInFile();
    Utils::separateur();

    return 0;
}
