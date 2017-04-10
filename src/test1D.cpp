#include <iostream>
#include <string>
#include <time.h>
#include "../include/LagrangeInterpolation1D.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    const string sep =  "\n********************************************************************************************************************\n\n";
    srand (time(NULL));

    // Get the size of the data points sequence from argument
    int size = (argc == 2) ? stoi(argv[1]) : Utils::randomValue(10,100) ;
    LagrangeInterpolation1D* interp = new LagrangeInterpolation1D(size);

    // Create the sequence of points uniformly
    interp->setPoints(Utils::createDataPoints(size));
    cout << "Sequence uniforme de " << interp->points().size() << " points d'interpolation " << endl;
    for (int i=0; i<int(interp->points().size()); ++i)
        cout << interp->points()[i] << " ";
    cout << endl << sep ;

    // Test de l'interpolation
    vector<double> testPoints;
    vector<double> realValue, estimate;
    double val = 0;
    testPoints.resize(int(Utils::randomValue(10,25)));
    cout << "Sequence de " << testPoints.size() << " points de test:" << endl;
    for (int i=0; i<int(testPoints.size()); i++)
    {
        testPoints[i] = Utils::randomValue(-1,1);
        cout << testPoints[i] << " ";
    }

    cout << endl << sep << "Calcul direct:" << endl;
    for (int i=0; i<int(testPoints.size()); i++)
    {
        val = interp->g(testPoints[i]);
        realValue.push_back(val);
        cout << val << " ";
    }
    interp->computeLiMinus1Fi(size);
    cout << endl << sep << "Calcul par interpolation:" << endl;
    for (int i=0; i<int(testPoints.size()); i++)
    {
        val = interp->lagrangeInterpolation_1D_iterative(testPoints[i],size);
        estimate.push_back(val);
        cout << val << " ";
    }
    cout << endl;
    cout << endl << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl << sep;

    return 0;
}
