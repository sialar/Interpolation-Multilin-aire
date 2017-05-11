#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include "../include/Interpolation.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    vector<int> nbPoints;
    int n = 9;
    nbPoints.push_back(n);

    Interpolation* interp = new Interpolation(nbPoints,1);
    interp->setDirPoints(0,nbPoints[0]);
    interp->displayInterpolationPoints();

    vector<double> x, y, z;
    double t;

    int function = -1;
    while (function !=0 && function != 1 && function!=2)
    {
        cout << "Choose the function to plot:" << endl;
        cout << " - 0: polynome de lagrange definit globalement" << endl;
        cout << " - 1: fonction affine par morceaux" << endl;
        cout << " - 2: polynome de lagrange definit par morceaux" << endl;

        cin >> function;
    }

    int point[2] = {3,4};
    int m = 100;
    double start, end;
    for (int i=0; i<m; i++)
    {
        start = -1 + 2.0*float(i)/float(m);
        end = -1 + 2.0*(float(i)+1.0)/float(m);
        t = Utils::randomValue(start, end);
        x.push_back(t);
        switch (function)
        {
            case 0:
                y.push_back(interp->lagrangeBasisFunction_1D(point[0],n,t,0));
                z.push_back(interp->lagrangeBasisFunction_1D(point[1],n,t,0));
                break;
            case 1:
                y.push_back(interp->piecewiseFunction_1D(point[0],t,0));
                z.push_back(interp->piecewiseFunction_1D(point[1],t,0));
                break;
            case 2:
                y.push_back(interp->piecewiseLagrangeBasisFunction_1D(point[0],n,t,0));
                z.push_back(interp->piecewiseLagrangeBasisFunction_1D(point[1],n,t,0));
                break;
            default:
                break;
        }
    }

    Utils::storeFunction(x,y,z);
    return 0;
}
