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
    MultiVariatePoint<int> methods(dim,-1);
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

string chooseCoreType(int argc, char* argv[], int argNum)
{
    if (argc > argNum)
    {
        if (!Functions::validCoreType(argv[argNum]))
        {
            cout << "Invalid core type!" << endl;
            exit(1);
        }

      return argv[argNum];
    }
    string core = "";
    while (!Functions::validCoreType(core))
    {
      cout << " - Choose the core type from { 'MOX', 'UOX', 'UOX-Gd' } : ";
      cin >> core;
    }
    Utils::separateur();
    return core;
}

string chooseCrossSection(int argc, char* argv[], int argNum)
{
    if (argc > argNum)
    {
        if (!Functions::validCrossSection(argv[argNum]))
        {
            cout << "Invalid cross section name!" << endl;
            exit(1);
        }

      return argv[argNum];
    }
    string cs = "";
    while (!Functions::validCrossSection(cs))
    {
      cout << " - Choose the cross section name from ";
      cout << "{ 'macro_totale0', 'macro_totale1', 'macro_absorption0', 'macro_absorption1', "
           << "'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', "
           << "'macro_nu*fission0', 'macro_nu*fission1', 'macro_fission0', 'macro_fission1' } : ";
      cin >> cs;
    }
    Utils::separateur();
    return cs;
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
  Utils::separateur();
  return maxIteration;
}

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    Utils::separateur();
    Functions::createFunctionsDataBase();
    int dimD = 5;
    int maxIteration = chooseMaxIteration(argc,argv,1);
    string core = chooseCoreType(argc,argv,2);
    string cs = chooseCrossSection(argc,argv,3);

    MultiVariatePoint<int> methods(5,0);// = chooseMethods(dimD);
    methods(0) = 2;
    MixedInterpolationPtr interp(new MixedInterpolation(dimD,core,cs,maxIteration,methods));
    //LagrangeInterpolationPtr interp(new LagrangeInterpolation(dimD,core,reactions,maxIteration));

    interp->readDataAndResults();
    interp->launchAIAlgo(false);
    interp->saveResults();
    interp->displayAll();

    return 0;
}
