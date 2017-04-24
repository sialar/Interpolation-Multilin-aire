#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/LagrangeInterpolationND.hpp"
#include "../include/Utils.hpp"

using namespace std;

vector<int> initData(int argc, char* argv[])
{
    int d = 0;
    vector<int> nbPoints;
    cout << endl;
    if (argc<2)
    {
        while (d<2 || d>10)
        {
            cout << "Choisir la dimension d (entre 2 et 10): ";
            cin >> d;
        }
        nbPoints.resize(d);
        for (int i=0; i<d; i++)
        {
            while (nbPoints[i]<1)
            {
                cout << "Choisir le nombre de points suivant la direction " << i << ": " ;
                cin >> nbPoints[i];
            }
        }
    }
    else
    {
        d = argc - 1;
        cout << "Dimension d: " << d << endl;
        nbPoints.resize(d);
        for (int i=0; i<d; i++)
        {
            nbPoints[i] = stoi(argv[i+1]);
            cout << "Nombre de points suivant la direction " << i << ": " << nbPoints[i] << endl;
        }
    }
    return nbPoints;
}
vector<double> createTestPointRandomly(int d)
{
    vector<double> testPoint;
    testPoint.resize(d);
    cout << "Point de test: (";
    for (int i=0; i<d; i++)
        testPoint[i] = Utils::randomValue(-1,1);
    for (int i=0; i<d-1; i++)
        cout << testPoint[i] << ",";
    cout << testPoint[d-1] << ")" << endl;
    return testPoint;
}

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    vector<int> nbPoints = initData(argc,argv);
    int d = nbPoints.size();

    LagrangeInterpolationND* interp = new LagrangeInterpolationND(nbPoints);
    Utils::separateur();

    for (int i=0; i<d; i++)
        interp->setDirPoints(i,Utils::createLejaSequence(nbPoints[i]));
    interp->showPoints();
    Utils::separateur();

    vector<double> testPoint = createTestPointRandomly(d);
    cout << Utils::gNd(testPoint) << endl;
    Utils::separateur();

    return 0;
}
