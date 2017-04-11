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
    const string sep =  "\n********************************************************************************************************************\n\n";

    // Get the size of the data points sequence from argument
    int size = (argc == 2) ? stoi(argv[1]) : Utils::randomValue(10,100) ;
    LagrangeInterpolation1D* interp = new LagrangeInterpolation1D(size);

    // Create the sequence of points uniformly
    interp->setPoints(Utils::createUniformSequence(size));
    cout << "Sequence uniforme de " << interp->points().size() << " points d'interpolation " << endl;
    for (int i=0; i<int(interp->points().size()); ++i)
        cout << interp->points()[i] << " ";
    cout << endl << sep ;

    // Test of the interpolation
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
    t1 = clock(); // start counting
    val = interp->lagrangeInterpolation_1D_iterative(Utils::randomValue(-1,1),size);
    t2 = clock(); // stop the count
    temps = (float)(t2-t1)/CLOCKS_PER_SEC; // compute the execution time
    cout << endl << "Le temps nécessaire pour interpoler la fonction g en un point choisi aléatoirement entre -1 et 1, est: "
                 << temps << ". (la valeur de g est connue en " << size << " points)" << endl << sep;

    cout << "L'erreur quadratique moyenne: " << Utils::squareError(realValue,estimate) << endl << sep;
    /*
    // Test Leja sequence
    vector<double> pp = Utils::createLejaSequence(10);
    for (size_t i=0; i<pp.size(); i++)
        cout << pp[i] << " ";
    cout << endl;
    */
    return 0;
}
