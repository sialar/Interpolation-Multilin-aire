#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include <limits>
#include <stdio.h>

#include "../include/Utils.hpp"
#include "../include/Functions.hpp"
#include "../include/MixedInterpolation.hpp"

using namespace std;

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
            return Functions::allCrossSectionType;
        for (int i=argNum; i<argc; i++)
        {
            if (!Functions::validCrossSections(argv[i]))
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
    if (nbReactions==12) return Functions::allCrossSectionType;
    for (int i=0; i<nbReactions; i++)
    {
        reaction = "";
        while (!Functions::validCrossSections(reaction))
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

void saveErrorsInFile(vector<vector<double>> results, string core, vector<string> reactions, vector<int> iterations)
{
    ofstream file(Utils::projectPath + "AI/data/" + core + "/RelativeErrors", ios::out );
    if(file)
    {
        file << results.size() << endl;
        for (int i=0; i<int(iterations.size()); i++)
            file << iterations[i] << " ";
        file << endl;
        file << reactions.size() << endl;
        for (int i=0; i<int(reactions.size()); i++)
            file << reactions[i] << " ";
        file << endl;
        for (int i=0; i<int(results.size()); i++)
        {
            for (int j=0; j<int(reactions.size()); j++)
                file << results[i][j] << " ";
            file << endl;
        }
        file.close();
    }
    else cerr << "Error while opening the file!" << endl;
}

vector<double> interpolate(int dimD,string core, vector<string> reactions, int nIter, MultiVariatePoint<int> methods)
{
    MixedInterpolationPtr interp(new MixedInterpolation(dimD,core,reactions,nIter,methods));

    interp->readDataAndResults();
    interp->launchAIAlgo(false);
    interp->saveResults();

    return interp->relativeErrors();
}

int main( int argc, char* argv[] )
{
    Functions::createFunctionsDataBase();
    vector<vector<double>> nIter_errors;
    int dimD = 5;
    vector<int> iterations(7);
    string core = chooseCoreType(argc,argv,1);
    vector<string> reactions = chooseReactionsType(argc,argv,2);
    MultiVariatePoint<int> methods(5,0,0);
    methods(0) = 2;

    iterations[0] = 50;
    iterations[1] = 100;
    iterations[2] = 200;
    iterations[3] = 300;
    iterations[4] = 450;
    iterations[5] = 600;
    iterations[6] = 800;

    for (int i=0; i<int(iterations.size()); i++)
    {
        nIter_errors.push_back(interpolate(dimD,core,reactions,iterations[i],methods));
        cout << "Computing relative errors with " << iterations[i] << " iterations done!" << endl;
    }

    saveErrorsInFile(nIter_errors, core, reactions, iterations);
    return 0;
}
