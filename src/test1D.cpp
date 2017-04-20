#include <iostream>
#include <string>
#include <time.h>
#include "../include/LagrangeInterpolation1D.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    // Lecture des données à partir les arguments
    int size = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(10,100);
    int nbTestPoints = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,100);
    LagrangeInterpolation1D* interp = new LagrangeInterpolation1D(size);

    // Evaluation de la fonction g (points choisis aléatoirement entre -1 et 1)
      vector<double> testPoints, realValues, estimate;
      testPoints.resize(nbTestPoints);
      for (int i=0; i<nbTestPoints; i++)
          testPoints[i] = Utils::randomValue(-1,1);

      Utils::displayTestPoints(testPoints,testPoints,1);
      Utils::separateur();
      realValues = Utils::displayGRealValues(testPoints,testPoints,1,true);
      Utils::separateur();


    /**************************************************************************/
    /**************** 1ere méthode : sequence uniforme ************************/
    /**************************************************************************/
    cout << "   - Methode d'interpolation utilisant une sequence uniforme de " << size << " points:" << endl;

    // Creation la sequence uniforme
    vector<double> uniformSequence = Utils::createUniformSequence(size);
    Utils::displayInterpolationPoints(uniformSequence,uniformSequence,1);


    // Test de l'interpolation en utilisant la sequence uniforme
    cout << endl << "   - Calcul par interpolation (evaluation de ĝ en " << nbTestPoints << " points de test)" << endl;
    interp->setPoints(uniformSequence);
    interp->computeAllAlphaI(size);
    for (int i=0; i<int(testPoints.size()); i++)
        estimate.push_back(interp->lagrangeInterpolation_1D_iterative(testPoints[i],size));
    Utils::displayApproximation(estimate,testPoints.size(),testPoints.size(),1,true);

    // Evaluation
    cout << endl << "   - Temps d'éxecution: " << interp->computeExecTimeOfOneApprox() << endl;
    cout << "   - Erreur quadratique moyenne: " << Utils::squareError(realValues,estimate) << endl;
    Utils::separateur();


    /**************************************************************************/
    /**************** 2eme méthode : sequence de Leja *************************/
    /**************************************************************************/
    cout << "   - Methode d'interpolation utilisant une sequence de " << size << " points de Leja:" << endl;
    estimate.clear();

    // Creation de la sequence de Leja
    vector<double> lejaSequence = Utils::createLejaSequence(size);
    Utils::displayInterpolationPoints(lejaSequence,lejaSequence,1);


    // Test de l'interpolation en utilisant la sequence de Leja
    cout << endl << "   - Calcul par interpolation (evaluation de ĝ en " << nbTestPoints << " points de test)" << endl;
    interp->setPoints(lejaSequence);
    interp->computeAllAlphaI(size);
    for (int i=0; i<int(testPoints.size()); i++)
        estimate.push_back(interp->lagrangeInterpolation_1D_iterative(testPoints[i],size));
    Utils::displayApproximation(estimate,testPoints.size(),testPoints.size(),1,true);

    // Sauvegarde des donnees pour affichage
    Utils::storeResult1D(testPoints,estimate,realValues);

    // Evaluation
    cout << endl << "   - Temps d'éxecution: " << interp->computeExecTimeOfOneApprox() << endl;
    cout << "   - Erreur quadratique moyenne: " << Utils::squareError(realValues,estimate) << endl;
    Utils::separateur();

    return 0;
}
