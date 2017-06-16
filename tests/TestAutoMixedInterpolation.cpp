#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/MixedInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    Utils::separateur();
    int dimD = Utils::chooseDimensionD(argc,argv,1);
    int dimN = Utils::chooseDimensionN(argc,argv,2);
    int nbTestPoints = Utils::chooseNbTestPoints(argc,argv,3);
    MultiVariatePoint<int> methods(dimD,0,0);
    int maxIteration = Utils::chooseMaxIteration(argc,argv,4);
    int f = Utils::chooseFunction(argc,argv,5);
    bool save = Utils::saveResults(argc,argv,6);

    Function interpFunc;
    if (f==1) interpFunc = Functions::autoPolynomialFunction;
    if (f==2) interpFunc = Functions::functionToPlot;
    if (f==3) interpFunc = Functions::sinOfNorm2;

    MixedInterpolationPtr interp(new MixedInterpolation(dimD,dimN,maxIteration,methods,interpFunc));
    interp->disableProgressDisplay();

    // Initialisation of test points
    interp->setRandomTestPoints(nbTestPoints);

    // Path creation
    double threshold = 1e-9;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    cout << " - The algorithm will be tested with different methods in each direction." << endl;
    cout << " - Each time, the algorithm will stop when the interpolation error becomes lower than"
         << " a threshold = " << threshold << endl;

    MultiVariatePoint<int> optimalMethods = interp->tryAllCases(threshold);
    interp->setMethods(optimalMethods);

    Utils::separateur();
    interp->tryWithDifferentMethods(optimalMethods, threshold);
    if (save) interp->savePathInFile(Utils::projectPath + "data/path.txt");
    Utils::separateur();

    return 0;
}
