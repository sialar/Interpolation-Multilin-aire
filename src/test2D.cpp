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
    const string sep =  "\n********************************************************************************************************************\n\n";

    // Get the size of the data points sequence from argument
    int sizeX = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(10,100);
    int sizeY = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,100);
    int version = (argc > 3) ? stoi(argv[3]) : 0;
    LagrangeInterpolation2D* interp = new LagrangeInterpolation2D(sizeX, sizeY, version);

    // Create the sequence of points uniformly
    interp->setPointsX(Utils::createUniformSequence(sizeX));
    interp->setPointsY(Utils::createUniformSequence(sizeY));

    cout << "Sequence uniforme des points d'interpolation: (" << interp->pointsX().size() <<
            " suivant l'axe des x et " <<  interp->pointsY().size() << " suivant l'axe des y)" << endl;
    for (int i=0; i<int(interp->pointsX().size()); ++i)
        cout << interp->pointsX()[i] << " ";
    cout << endl;
    for (int j=0; j<int(interp->pointsY().size()); ++j)
        cout << interp->pointsY()[j] << " ";
    cout << endl << sep;

    // Test of the interpolation
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
/*
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

    cout << endl << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl;
*/

    cout << sep << "Calcul par interpolation bilinéaire:" << endl;
    interp->computeValues(sizeX,sizeY);

    indice2D lastIndexInPath = interp->m_path[interp->m_path.size()-1];
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
    temps = (float)(t2-t1)/CLOCKS_PER_SEC; // compute the execution time
    cout << endl << "Le temps nécessaire pour interpoler la fonction g en un point 2d dont les composantes sont choisies aléatoirement entre -1 et 1, est: "
                 << temps << ". (la valeur de g est connue en " << sizeX << "x" << sizeY << "=" << sizeX*sizeY << " points)" << endl << sep;

    cout << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl << sep;
    return 0;
}
