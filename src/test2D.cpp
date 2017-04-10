#include <iostream>
#include <string>
#include <time.h>
#include "../include/LagrangeInterpolation2D.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    const string sep =  "\n********************************************************************************************************************\n\n";
    srand (time(NULL));

    // Get the size of the data points sequence from argument
    int sizeX = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(10,100) ;
    int sizeY = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,100) ;
    LagrangeInterpolation2D* interp = new LagrangeInterpolation2D(sizeX, sizeY);

    // Create the sequence of points uniformly
    interp->setPointsX(Utils::createDataPoints(sizeX));
    interp->setPointsY(Utils::createDataPoints(sizeY));

    cout << "Sequence uniforme des points d'interpolation: (" << interp->pointsX().size() <<
            " suivant l'axe des x et " <<  interp->pointsY().size() << " suivant l'axe des y)" << endl;
    for (int i=0; i<int(interp->pointsX().size()); ++i)
        cout << interp->pointsX()[i] << " ";
    cout << endl;
    for (int j=0; j<int(interp->pointsY().size()); ++j)
        cout << interp->pointsY()[j] << " ";
    cout << endl << sep;

    // Test de l'interpolation
    vector<double> testPointsX, testPointsY;
    vector<double> realValue, estimate;
    double val = 0;
    testPointsX.resize(int(Utils::randomValue(10,15)));
    testPointsY.resize(int(Utils::randomValue(10,15)));
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
    cout << endl << sep << "Calcul direct:" << endl;
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

    cout << sep << "Calcul par une simple interpolation bilinéaire:" << endl;
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->lagrangeInterpolation_2D_simple(testPointsX[i],testPointsY[j]);
            estimate.push_back(val);
            cout << val << " ";
        }
        cout << endl;
    }

    cout << endl << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl << sep;
    return 0;
}
