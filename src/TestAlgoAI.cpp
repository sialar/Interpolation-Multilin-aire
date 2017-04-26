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
    while (d<2 || d>10)
    {
        cout << " - Choose the space dimension (between 2 and 10) : ";
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

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    int nbTestPoints = 0;
    vector<MultiVariatePoint<double>> testPoints;
    vector<double> realValues, estimate;

    // Initialisation of interpolation points
    vector<int> nbPoints = initData(argc,argv);
    Interpolation* interp = new Interpolation(nbPoints);
    Utils::separateur();
    for (size_t i=0; i<nbPoints.size(); i++)
        interp->setDirPoints(i,Utils::createLejaSequence(nbPoints[i]));
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
    cout << " - Sequence of " << nbTestPoints << " random test point : " ;
    Utils::displayPoints(testPoints);
    Utils::separateur();

    // Path creation
    int nbIteration = 0;
    int maxIteration = 1;
    for (size_t i=0; i<nbPoints.size(); i++)
        maxIteration *= interp->points()[i].size();
    while (nbIteration<1)
    {
        cout << " - Choose the number of iteration ( < " << maxIteration << " ) in AI algo : " ;
        cin >> nbIteration;
    }

    interp->testPathBuilt(nbIteration);
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
        estimate.push_back(interp->lagrangeInterpolation_ND_iterative(p));
    Utils::displayPoints(estimate);

    // Evaluation
    cout << endl << " - Mean squared error = " << Utils::squareError(realValues,estimate) << endl;
    Utils::separateur();

    return 0;
}
