#include <iostream>
#include <string>
#include <time.h>
#include "../include/LagrangeInterpolation2D.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    float temps;
    clock_t t1, t2;


    // Get the size of the data points sequence from argument
    int sizeX = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(10,100);
    int sizeY = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,100);
    int version = (argc > 3) ? stoi(argv[3]) : 0;
    LagrangeInterpolation2D* interp = new LagrangeInterpolation2D(sizeX, sizeY, version);


    // Evaluation de la fonction g (points choisis aléatoirement entre -1 et 1)
    vector<double> testPointsX, testPointsY;
    vector<double> realValue, estimate;
    double val = 0;
    testPointsX.resize(int(Utils::randomValue(2,5)));
    testPointsY.resize(int(Utils::randomValue(2,5)));
    cout << "Sequence des points de test: (" << testPointsX.size() <<
            " suivant l'axe des x et " << testPointsY.size() << " suivant l'axe des y)" << endl;
    cout << "-  Suivant l'axe des x: ";
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        testPointsX[i] = Utils::randomValue(-1,1);
        cout << testPointsX[i] << " ";
    }
    cout << endl << "-  Suivant l'axe des y: ";
    for (int j=0; j<int(testPointsY.size()); j++)
    {
        testPointsY[j] = Utils::randomValue(-1,1);
        cout << testPointsY[j] << " ";
    }
    cout << endl << endl << "Calcul direct:" << endl;
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->g(testPointsX[i],testPointsY[j]);
            realValue.push_back(val);
            cout << val << " ";
        }
        cout << endl;
    }

    Utils::separateur();
    interp->showPath();

    /**************************************************************************/
    /**************** 1ere méthode : sequence uniforme ************************/
    /**************************************************************************/
    Utils::separateur();
    cout << "   - 1ere methode:" << endl;
    // Creation de la sequence uniforme
    interp->setPointsX(Utils::createUniformSequence(sizeX));
    interp->setPointsY(Utils::createUniformSequence(sizeY));
    cout << endl << "Sequence uniforme (" << interp->pointsX().size() <<
            "x" <<  interp->pointsY().size() << " points)" << endl;
    for (int i=0; i<int(interp->pointsX().size()); ++i)
        cout << interp->pointsX()[i] << " ";
    cout << endl;
    for (int j=0; j<int(interp->pointsY().size()); ++j)
        cout << interp->pointsY()[j] << " ";
    cout << endl;
    // Test de l'interpolation en utilisant la sequence uniforme
    cout << endl << "Calcul par interpolation bilinéaire:" << endl;
    interp->computeAllAlphaNu(sizeX,sizeY);
    indice2D lastIndexInPath = interp->path()[interp->path().size()-1];
    estimate.clear();
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j],lastIndexInPath[0],lastIndexInPath[1]);
            estimate.push_back(val);
            cout << val << " ";
        }
        cout << endl;
    }
    t1 = clock(); // start counting
    val = interp->lagrangeInterpolation_2D_iterative(Utils::randomValue(-1,1),Utils::randomValue(-1,1),lastIndexInPath[0],lastIndexInPath[1]);
    t2 = clock(); // stop the count
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << endl << "Temps d'éxecution: " << temps << endl;
    cout << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl;


    /**************************************************************************/
    /**************** 2eme méthode : sequence de Leja *************************/
    /**************************************************************************/

    Utils::separateur();
    estimate.clear();
    cout << "   - 2eme methode:" << endl;
    // Creation de la sequence de Leja
    interp->setPointsX(Utils::createLejaSequence(sizeX));
    interp->setPointsY(Utils::createLejaSequence(sizeY));
    cout << endl << "Sequence de Leja (" << interp->pointsX().size() <<
            "x" <<  interp->pointsY().size() << " points)" << endl;
    for (int i=0; i<int(interp->pointsX().size()); ++i)
        cout << interp->pointsX()[i] << " ";
    cout << endl;
    for (int j=0; j<int(interp->pointsY().size()); ++j)
        cout << interp->pointsY()[j] << " ";
    cout << endl;
    // Test de l'interpolation en utilisant la sequence de Leja
    cout << endl << "Calcul par interpolation bilinéaire:" << endl;
    interp->computeAllAlphaNu(sizeX,sizeY);
    lastIndexInPath = interp->path()[interp->path().size()-1];
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j],lastIndexInPath[0],lastIndexInPath[1]);
            estimate.push_back(val);
            cout << val << " ";
        }
        cout << endl;
    }
    t1 = clock(); // start counting
    val = interp->lagrangeInterpolation_2D_iterative(Utils::randomValue(-1,1),Utils::randomValue(-1,1),lastIndexInPath[0],lastIndexInPath[1]);
    t2 = clock(); // stop the count
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << endl << "Temps d'éxecution: " << temps << endl;
    cout << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl;

    /**************************************************************************/
    /************** 3eme méthode : sequence de Leja et algo AI ****************/
    /**************************************************************************/

    Utils::separateur();
    estimate.clear();
    interp->buildPathWithAIAlgo(sizeX,sizeY,true);
    cout << "Chemin obtenu avec l'algo AI:" << endl;
    interp->showPath();
    cout << endl << "Calcul par interpolation bilinéaire en utilisant l'algorithme AI:" << endl;
    //interp->computeAllAlphaNu(sizeX,sizeY);
    lastIndexInPath = interp->path()[interp->path().size()-1];
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j],lastIndexInPath[0],lastIndexInPath[1]);
            estimate.push_back(val);
            cout << val << " ";
        }
        cout << endl;
    }
    t1 = clock(); // start counting
    val = interp->lagrangeInterpolation_2D_iterative(Utils::randomValue(-1,1),Utils::randomValue(-1,1),lastIndexInPath[0],lastIndexInPath[1]);
    t2 = clock(); // stop the count
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << endl << "Temps d'éxecution: " << temps << endl;
    cout << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl;
    return 0;
}
