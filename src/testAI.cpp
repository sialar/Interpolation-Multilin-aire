#include <iostream>
#include <string>
#include <time.h>
#include "../include/LagrangeInterpolation2D.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    // Lecture des données à partir les arguments
    int n = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(10,100);
    int m = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,100);

    int nbTestPointsX = (argc > 3) ? stoi(argv[3]) : Utils::randomValue(2,5);
    int nbTestPointsY = (argc > 4) ? stoi(argv[4]) : Utils::randomValue(2,5);
    int nbIteration = (argc > 5) ? stoi(argv[5]) : n*m;
    LagrangeInterpolation2D* interp = new LagrangeInterpolation2D(n, m, -1);

    // Evaluation de la fonction g (points choisis aléatoirement entre -1 et 1)
    vector<double> testPointsX, testPointsY, realValues, estimate;
    Utils::separateur();
    testPointsX.resize(nbTestPointsX);
    for (int i=0; i<nbTestPointsX; i++)
        testPointsX[i] = Utils::randomValue(-1,1);
    testPointsY.resize(nbTestPointsY);
    for (int i=0; i<nbTestPointsY; i++)
        testPointsY[i] = Utils::randomValue(-1,1);

    Utils::displayTestPoints(testPointsX,testPointsY,2);
    Utils::separateur();
    realValues = Utils::displayGRealValues(testPointsX,testPointsY,2,true);
    Utils::separateur();

    interp->setPointsX(Utils::createLejaSequence(n));
    interp->setPointsY(Utils::createLejaSequence(m));

    interp->initAlphaTab(numeric_limits<double>::max());
    //interp->buildPathWithAIAlgo(nbIteration);
    interp->testPathBuilt(nbIteration);

    if (n * m <= 100)
    {
        cout << "   - Chemin obtenu avec l'algo AI:" << endl;
        interp->showPath();
        Utils::separateur();
    }
    interp->savePathInFile();

    // Approximer g
    for (int i=0; i<int(testPointsX.size()); i++)
        for (int j=0; j<int(testPointsY.size()); j++)
          estimate.push_back(interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j]));
    Utils::displayApproximation(estimate,testPointsX.size(),testPointsY.size(),2,true);

    // Evaluation
    cout << endl << "   - Temps d'éxecution: " << interp->computeExecTimeOfOneApprox() << endl;
    cout << "   - Erreur quadratique moyenne: " << Utils::squareError(realValues,estimate) << endl;
    Utils::separateur();

    return 0;
}
