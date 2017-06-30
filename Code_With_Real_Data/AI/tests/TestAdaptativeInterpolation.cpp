#include <python3.5/Python.h>
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

vector<string> chooseReactionsType(int argc, char* argv[], int argNum)
{
    vector<string> reactions;
    if (argc > argNum)
    {
        if (string(argv[argNum]).compare("ALL")==0)
            return Functions::allReactionTypes;
        for (int i=argNum; i<argc; i++)
        {
            if (!Functions::validReactionType(argv[i]))
            {
                cout << "Invalid reaction type!" << endl;
                exit(1);
            }
            reactions.push_back(argv[i]);
        }
        return reactions;
    }
    int nbReactions = -1;
    string reaction = "";
    while (nbReactions < 0 || nbReactions > 12)
    {
      cout << " - Choose the number of reaction types : ";
      cin >> nbReactions;
    }
    if (nbReactions==12) return Functions::allReactionTypes;
    for (int i=0; i<nbReactions; i++)
    {
        reaction = "";
        while (!Functions::validReactionType(reaction))
        {
          cout << " - Choose the reaction type " << i+1 << " from { 'macro_totale0', 'macro_totale1', "
               << "'macro_absorption0', 'macro_absorption1', 'macro_scattering000', "
               << "'macro_scattering001', 'macro_scattering010', 'macro_scattering011', "
               << "'macro_nu*fission0', 'macro_nu*fission1', 'macro_fission0', 'macro_fission1', 'ALL' }\n - ";
          cin >> reaction;
        }
        reactions.push_back(reaction);
    }
    Utils::separateur();
    return reactions;
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
    vector<string> reactions = chooseReactionsType(argc,argv,3);

    //MultiVariatePoint<int> methods = chooseMethods(dimD);
    //MixedInterpolationPtr interp(new MixedInterpolation(dimD,reactions.size(),maxIteration,methods));

    LagrangeInterpolationPtr interp(new LagrangeInterpolation(dimD,reactions.size(),maxIteration));
    interp->setFunc(core,reactions);

    interp->readEDFTestPointsFromFile();
    interp->displayRealDomain();
    Utils::separateur();

    // Path creation
    double threshold = 1e-9;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = " \
         << threshold << endl;
    interp->buildPathWithAIAlgo(threshold, false);


    // Computing real values, and approximation of function g at test points
    Utils::separateur();
    vector<vector<double>> realValues, estimate;
    for (int i=0; i<10; i++)
    {
        realValues.push_back(interp->func(interp->testPoints()[i]));
        estimate.push_back(interp->interpolation(interp->testPoints()[i],interp->path().size()));
    }

    // Evaluation
    double relativeError = Utils::relativeInterpolationError(realValues,estimate);
    double mseError = Utils::mseInterpolationError(realValues,estimate);
    cout << " - Relative Interpolation error (pcm) = " << relativeError << endl;
    cout << " - MSE Interpolation error (pcm) = " << mseError << endl;
    cout << " - Number of evaluation = " << interp->nbEvals() << endl;
    cout << " - Total Time = " << interp->totalTime() << endl;
    cout << " - AI Run Time = " << interp->runTime() << endl;

    return 0;
}
