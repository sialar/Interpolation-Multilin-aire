#include <iostream>
#include <string>
#include <time.h>
#include "../include/LagrangeInterpolation1D.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    float temps;
    clock_t t1, t2;


    // Get the size of the data points sequence from argument
    int size = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(10,100);
    int nbTestPoints = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,100);
    LagrangeInterpolation1D* interp = new LagrangeInterpolation1D(size);


    // Evaluation de la fonction g (points choisis aléatoirement entre -1 et 1)
    vector<double> testPoints;
    vector<double> realValue, estimate;
    double val = 0;
    testPoints.resize(nbTestPoints);
    cout << endl << "Sequence de " << testPoints.size() << " points de test:" << endl;
    for (int i=0; i<int(testPoints.size()); i++)
    {
        testPoints[i] = Utils::randomValue(-1,1);
        cout << testPoints[i] << " ";
    }
    cout << endl << endl << "Calcul direct:" << endl;
    for (int i=0; i<int(testPoints.size()); i++)
    {
        val = interp->g(testPoints[i]);
        realValue.push_back(val);
        cout << val << " ";
    }
    cout << endl;


    /**************************************************************************/
    /**************** 1ere méthode : sequence uniforme ************************/
    /**************************************************************************/

    Utils::separateur();
    cout << "  ---> 1ere methode:" << endl;
    // Creation la sequence uniforme
    vector<double> uniformSequence = Utils::createUniformSequence(size);
    cout << endl << "Sequence uniforme (" << size << " points)" << endl;
    for (int i=0; i<size; ++i)
        cout << uniformSequence[i] << " ";
    // Test de l'interpolation en utilisant la sequence uniforme
    cout << endl << endl << "Calcul par interpolation:" << endl;
    interp->setPoints(uniformSequence);
    interp->computeAllAlphaI(size);
    for (int i=0; i<int(testPoints.size()); i++)
    {
        val = interp->lagrangeInterpolation_1D_iterative(testPoints[i],size);
        estimate.push_back(val);
        cout << val << " ";
    }
    cout << endl;
    t1 = clock(); // start counting
    val = interp->lagrangeInterpolation_1D_iterative(Utils::randomValue(-1,1),size);
    t2 = clock(); // stop the count
    temps = (float)(t2-t1)/CLOCKS_PER_SEC; // compute the execution time
    cout << endl << "Temps d'éxecution: " << temps << endl;
    cout << endl << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl;


    /**************************************************************************/
    /**************** 2eme méthode : sequence de Leja *************************/
    /**************************************************************************/

    Utils::separateur();
    cout << "  ---> 2eme methode:" << endl;
    estimate.clear();
    // Creation de la sequence de Leja
    vector<double> lejaSequence = Utils::createLejaSequence(size);
    cout << endl << "Sequence de Leja (" << size << " points) " << endl;
    for (int i=0; i<size; ++i)
        cout << lejaSequence[i] << " ";
    // Test de l'interpolation en utilisant la sequence de Leja
    cout << endl << endl << "Calcul par interpolation:" << endl;
    interp->setPoints(lejaSequence);
    interp->computeAllAlphaI(size);
    for (int i=0; i<int(testPoints.size()); i++)
    {
        val = interp->lagrangeInterpolation_1D_iterative(testPoints[i],size);
        estimate.push_back(val);
        cout << val << " ";
    }
    cout << endl;
    t1 = clock(); // start counting
    val = interp->lagrangeInterpolation_1D_iterative(Utils::randomValue(-1,1),size);
    t2 = clock(); // stop the count
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << endl << "Temps d'éxecution: " << temps << endl;
    cout << endl << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl;


    return 0;
}
