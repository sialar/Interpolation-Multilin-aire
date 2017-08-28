#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/MixedInterpolation.hpp"
#include "../include/LagrangeInterpolation.hpp"
#include "../include/AnalyticalFunctions.hpp"
#include "../include/Utils.hpp"

using namespace std;

int chooseDimensionD(int argc, char* argv[], int argNum)
{
    if (argc > argNum) return stoi(argv[argNum]);
    int dim = -1;
    while (dim < 0)
    {
        cout << " - Choose the space dimension d: ";
        cin >> dim;
    }
    return dim;
}

int chooseDimensionN(int argc, char* argv[], int argNum)
{
    if (argc > argNum) return stoi(argv[argNum]);
    int dim = -1;
    while (dim < 0)
    {
        cout << " - Choose the space dimension n: ";
        cin >> dim;
    }
    return dim;
}

int chooseNbTestPoints(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  int nbTestPoints = -1;
  while (nbTestPoints < 0)
  {
    cout << " - Choose the number ot test points : ";
    cin >> nbTestPoints;
  }
  return nbTestPoints;
}

int chooseMaxIteration(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  int maxIteration = -1;
  while (maxIteration < 0)
  {
    cout << " - Choose the maximum number of iteration : ";
    cin >> maxIteration;
  }
  return maxIteration;
}

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
  	int d = chooseDimensionD(argc, argv, 1);
  	int n = chooseDimensionD(argc, argv, 2);
    AnalyticalFunctionsPtr f = make_shared<AnalyticalFunctions>(d,n);

  	int nbIter = chooseMaxIteration(argc, argv, 3);
  	int nbTestPts = chooseNbTestPoints(argc, argv, 4);
    MultiVariatePoint<int> methods(5,0,0);// = chooseMethods(d);

    MixedInterpolationPtr interp(new MixedInterpolation(f,nbIter,methods));
	  interp->setRandomTestPoints(nbTestPts);

    interp->launchAIAlgo(false);
    interp->computeAIApproximationResults();

	   interp->displayAll();

    return 0;
}
