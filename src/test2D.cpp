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
    int version = (argc > 5) ? stoi(argv[5]) : 1;
    int nbTestPointsX = (argc > 3) ? stoi(argv[3]) : Utils::randomValue(2,5);
    int nbTestPointsY = (argc > 4) ? stoi(argv[4]) : Utils::randomValue(2,5);
    LagrangeInterpolation2D* interp = new LagrangeInterpolation2D(n, m, version);
    bool debug1 =  nbTestPointsX*nbTestPointsY<101;
    bool debug2 =  interp->path().size()<101;

    // Evaluation de la fonction g (points choisis aléatoirement entre -1 et 1)
    vector<double> testPointsX, testPointsY, realValues, estimate;
    Utils::separateur();
    testPointsX.resize(nbTestPointsX);
    for (int i=0; i<nbTestPointsX; i++)
        testPointsX[i] = Utils::randomValue(-1,1);
    testPointsY.resize(nbTestPointsY);
    for (int i=0; i<nbTestPointsY; i++)
        testPointsY[i] = Utils::randomValue(-1,1);


    vector<vector<double>> testPointsSeq;
    testPointsSeq.push_back(testPointsX);
    testPointsSeq.push_back(testPointsY);
    cout << "   - Sequence de points de test:" << endl;
    Utils::displayPoints(testPointsSeq,2);
    Utils::separateur();

    realValues = Utils::displayGRealValues(testPointsX,testPointsY,2,debug1);
    Utils::separateur();
    if (debug2)
    {
        interp->showPath();
        Utils::separateur();
    }

    /**************************************************************************/
    /**************** 1ere méthode : sequence uniforme ************************/
    /**************************************************************************/
    estimate.clear();
    // Creation de la sequence uniforme
    interp->setPointsX(Utils::createUniformSequence(n));
    interp->setPointsY(Utils::createUniformSequence(m));
    cout << "   - Methode d'interpolation utilisant une sequence uniforme (" << n <<
            " points suivant la direction x et " << m << " points suivant la direction y)" << endl;
    vector<vector<double>> uniformSequences;
    uniformSequences.push_back(interp->pointsX());
    uniformSequences.push_back(interp->pointsY());
    Utils::displayPoints(uniformSequences,2);


    // Test de l'interpolation en utilisant la sequence uniforme
    cout << endl << "   - Calcul par interpolation bilinéaire: (evaluation de ĝ en " << nbTestPointsX*nbTestPointsY << " points de test)" << endl;
    interp->computeAllAlphaNu();
    for (int i=0; i<int(testPointsX.size()); i++)
        for (int j=0; j<int(testPointsY.size()); j++)
          estimate.push_back(interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j]));
    Utils::displayApproximation(estimate,testPointsX.size(),testPointsY.size(),2,debug1);

    // Evaluation
    cout << endl << "   - Temps d'éxecution: " << interp->computeExecTimeOfOneApprox() << endl;
    cout << "   - Erreur quadratique moyenne: " << Utils::squareError(realValues,estimate) << endl;
    Utils::separateur();

    /**************************************************************************/
    /**************** 2eme méthode : sequence de Leja *************************/
    /**************************************************************************/
    estimate.clear();
    // Creation de la sequence de Leja
    interp->setPointsX(Utils::createLejaSequence(n));
    interp->setPointsY(Utils::createLejaSequence(m));
    Utils::store2DLejaSequenceInFile(interp->pointsX(),interp->pointsY());
    cout << "   - Methode d'interpolation utilisant la sequence de Leja (" << interp->pointsX().size() <<
            " points suivant la direction x et " << interp->pointsY().size() << " points suivant la direction y)" << endl;
    vector<vector<double>> lejaSequences;
    lejaSequences.push_back(interp->pointsX());
    lejaSequences.push_back(interp->pointsY());
    Utils::displayPoints(lejaSequences,2);

    // Test de l'interpolation en utilisant la sequence de Leja
    cout << endl << "   - Calcul par interpolation bilinéaire: (evaluation de ĝ en " << nbTestPointsX*nbTestPointsY << " points de test)" << endl;
    interp->computeAllAlphaNu();
    for (int i=0; i<int(testPointsX.size()); i++)
        for (int j=0; j<int(testPointsY.size()); j++)
          estimate.push_back(interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j]));
    Utils::displayApproximation(estimate,testPointsX.size(),testPointsY.size(),2,debug1);

    // Evaluation
    cout << endl << "   - Temps d'éxecution: " << interp->computeExecTimeOfOneApprox() << endl;
    cout << "   - Erreur quadratique moyenne: " << Utils::squareError(realValues,estimate) << endl;
    Utils::separateur();

    /**************************************************************************/
    /************** 3eme méthode : sequence de Leja et algo AI ****************/
    /**************************************************************************/

    return 0;
}
