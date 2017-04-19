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
    int nbTestPointsX = (argc > 3) ? stoi(argv[3]) : Utils::randomValue(2,5);
    int nbTestPointsY = (argc > 4) ? stoi(argv[4]) : Utils::randomValue(2,5);
    int version = (argc > 5) ? stoi(argv[5]) : 1;
    LagrangeInterpolation2D* interp = new LagrangeInterpolation2D(sizeX, sizeY, version);
    bool debug =  (interp->path().size()<101) && (nbTestPointsX*nbTestPointsY<101);

    // Evaluation de la fonction g (points choisis aléatoirement entre -1 et 1)
    vector<double> testPointsX, testPointsY, realValue, estimate;
    double val = 0;
    testPointsX.resize(nbTestPointsX);
    testPointsY.resize(nbTestPointsY);
    Utils::separateur();
    cout << "   - Sequence des points de test: (" << testPointsX.size() <<
            " suivant l'axe des x et " << testPointsY.size() << " suivant l'axe des y):" << endl;
    cout << "         + Suivant la direction x: ";
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        testPointsX[i] = Utils::randomValue(-1,1);
        cout << testPointsX[i] << " ";
    }
    cout << endl << "         + Suivant la direction y: ";
    for (int j=0; j<int(testPointsY.size()); j++)
    {
        testPointsY[j] = Utils::randomValue(-1,1);
        cout << testPointsY[j] << " ";
    }
    cout << endl ;
    Utils::separateur();
    cout << "   - Calcul direct (evaluation de g en " << nbTestPointsX*nbTestPointsY << " points):" << endl;
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->g(testPointsX[i],testPointsY[j]);
            realValue.push_back(val);
            if (debug) cout << val << " ";
        }
        if (debug) cout << endl;
    }

    Utils::separateur();
    if (debug)
    {
        interp->showPath();
        Utils::separateur();
    }

    /**************************************************************************/
    /**************** 1ere méthode : sequence uniforme ************************/
    /**************************************************************************/

    estimate.clear();
    // Creation de la sequence uniforme
    interp->setPointsX(Utils::createUniformSequence(sizeX));
    interp->setPointsY(Utils::createUniformSequence(sizeY));
    cout << "   - Methode d'interpolation utilisant une sequence uniforme (" << interp->pointsX().size() <<
            " points suivant la direction x et " << interp->pointsY().size() << " points suivant la direction y)" << endl;
    cout << "         + Suivant la direction x: ";
    for (int i=0; i<int(interp->pointsX().size()); ++i)
        cout << interp->pointsX()[i] << " ";
    cout << endl << "         + Suivant la direction y: ";
    for (int j=0; j<int(interp->pointsY().size()); ++j)
        cout << interp->pointsY()[j] << " ";
    cout << endl;

    // Test de l'interpolation en utilisant la sequence uniforme
    cout << endl << "   - Calcul par interpolation bilinéaire: (evaluation de ĝ en " << nbTestPointsX*nbTestPointsY << " points de test)" << endl;
    interp->computeAllAlphaNu();
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j]);
            estimate.push_back(val);
            if (debug) cout << val << " ";
        }
        if (debug) cout << endl;
    }
    t1 = clock(); // start counting
    val = interp->lagrangeInterpolation_2D_iterative(Utils::randomValue(-1,1),Utils::randomValue(-1,1));
    t2 = clock(); // stop the count
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << endl << "   - Temps d'éxecution: " << temps << endl;
    cout << "   - Erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl;
    Utils::separateur();

    /**************************************************************************/
    /**************** 2eme méthode : sequence de Leja *************************/
    /**************************************************************************/
    estimate.clear();
    // Creation de la sequence de Leja
    interp->setPointsX(Utils::createLejaSequence(sizeX));
    interp->setPointsY(Utils::createLejaSequence(sizeY));
    Utils::store2DLejaSequenceInFile(interp->pointsX(),interp->pointsY());

    cout << "   - Methode d'interpolation utilisant la sequence de Leja (" << interp->pointsX().size() <<
            " points suivant la direction x et " << interp->pointsY().size() << " points suivant la direction y)" << endl;
    cout << "         + Suivant la direction x: ";
    for (int i=0; i<int(interp->pointsX().size()); ++i)
        cout << interp->pointsX()[i] << " ";
    cout << endl << "         + Suivant la direction y: ";
    for (int j=0; j<int(interp->pointsY().size()); ++j)
        cout << interp->pointsY()[j] << " ";
    cout << endl;

    // Test de l'interpolation en utilisant la sequence de Leja
    cout << endl << "   - Calcul par interpolation bilinéaire: (evaluation de ĝ en " << nbTestPointsX*nbTestPointsY << " points de test)" << endl;
    interp->computeAllAlphaNu();
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j]);
            estimate.push_back(val);
            if (debug) cout << val << " ";
        }
        if (debug) cout << endl;
    }
    Utils::storeResult2D(testPointsX,testPointsY,estimate,realValue);
    t1 = clock();
    val = interp->lagrangeInterpolation_2D_iterative(Utils::randomValue(-1,1),Utils::randomValue(-1,1));
    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << endl << "   - Temps d'éxecution: " << temps << endl;
    cout << "   - Erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl;
    Utils::separateur();

    /**************************************************************************/
    /************** 3eme méthode : sequence de Leja et algo AI ****************/
    /**************************************************************************/

    return 0;
}
