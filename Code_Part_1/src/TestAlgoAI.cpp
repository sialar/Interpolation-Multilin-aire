#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/Interpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;

vector<int> initData(int argc, char* argv[])
{
    int d = 0;
    vector<int> nbPoints;
    Utils::separateur();
    while (d<1 || d>20)
    {
        cout << " - Choose the space dimension (between 1 and 20) : ";
        cin >> d;
    }
    nbPoints.resize(d);
    for (int i=0; i<d; i++)
        while (nbPoints[i]<1)
        {
            cout << " - Choose the number of points in direction " << i << " : " ;
            cin >> nbPoints[i];
        }
    return nbPoints;
}

int chooseMethod()
{
    Utils::separateur();
    int method = -1;
    while (method!=1 && method!=0 && method!=2)
    {
        cout << " - Choose the method of interpolation: " << endl;
        cout << "\t - 0: Using lagrange polynomial functions and Leja points " << endl;
        cout << "\t - 1: Using fine-tune function functions and middle points: " << endl;
        cout << "\t - 2: Using piecewise lagrange polynomial functions and middle points: " << endl;
        cin >> method;
    }
    return method;
}

int main( int argc, char* argv[] )
{

    srand (time(NULL));
    int nbTestPoints = 0;
    vector<MultiVariatePoint<double>> testPoints;
    vector<double> realValues, estimate;

    // Initialisation of interpolation points
    vector<int> nbPoints = initData(argc,argv);
    Interpolation* interp = new Interpolation(nbPoints,chooseMethod());
    Utils::separateur();
    for (size_t i=0; i<nbPoints.size(); i++)
        interp->setDirPoints(i,nbPoints[i]);
    interp->displayInterpolationPoints();
    Utils::separateur();

    // Initialisation of test points
    while (nbTestPoints<1)
    {
        cout << " - Choose the number of test points : " ;
        cin >> nbTestPoints;
    }
    testPoints.resize(nbTestPoints);
    for (int j=0; j<nbTestPoints; j++)
        testPoints[j] = Utils::createRandomMultiVariatePoint(nbPoints.size());
    interp->setTestPoints(testPoints);
    cout << " - Sequence of " << nbTestPoints << " random test points : " ;
    Utils::displayPoints(testPoints);
    Utils::separateur();

    // Path creation
    long int maxIteration = 1;
    double threshold = 1e-4;
    for (size_t i=0; i<nbPoints.size(); i++)
        maxIteration *= interp->points()[i].size();

    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will stop when the interpolation error becomes lower than a threshold = "
         << threshold << endl << endl;
    interp->testPathBuilt(maxIteration, threshold, false);
    interp->savePathInFile();
    Utils::separateur();

    // Computing real values at test points
    cout << " - Real values of function g evaluated at test points :" << endl;
    for (MultiVariatePoint<double> p : testPoints)
        realValues.push_back(Utils::gNd(p));
    Utils::displayPoints(realValues);

    // Approximating g at test points
    cout << endl << " - Approximation of function g at test points : " << endl;
    for (MultiVariatePoint<double> p : testPoints)
        estimate.push_back(interp->interpolation_ND(p));
    Utils::displayPoints(estimate);

    // Evaluation
    cout << endl << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;
    Utils::separateur();

  //  cout << interp->lagrangeBasisFunction_1D(5,5,-0.577316 ,1) << endl;
  //  cout << interp->piecewiseFunction_1D(5,-0.577316,1) << endl;

    return 0;
}
