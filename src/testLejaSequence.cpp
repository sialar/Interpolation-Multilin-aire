#include <iostream>
#include <string>
#include <time.h>
#include "../include/LagrangeInterpolation2D.hpp"
#include "../include/Utils.hpp"
#include "../include/IndiceND.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    int sizeX = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(10,100);
    int sizeY = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,100);
    int sizeZ = (argc > 3) ? stoi(argv[3]) : Utils::randomValue(10,100);

    vector<double> pointsX = Utils::createLejaSequence(sizeX);
    vector<double> pointsY = Utils::createLejaSequence(sizeY);
    vector<double> pointsZ = Utils::createLejaSequence(sizeZ);

    if (argc <= 3)
        Utils::store2DLejaSequenceInFile(pointsX,pointsY);
    else
        Utils::store3DLejaSequenceInFile(pointsX,pointsY,pointsZ);

    return 0;
}
