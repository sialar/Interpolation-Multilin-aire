#include <iostream>
#include <string>
#include <time.h>
#include "../include/Utils.hpp"
#include "../include/MultiVariatePoint.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));

    int dim = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(2,3);
    if ((dim!=2) && (dim!=3))
    {
        cout << "Invalid argument 1: The space dimension must be in 2 or 3." << endl;
        exit(1);
    }

    int size = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,20);
    if ((size>1000) || (size<10))
    {
        cout << "Invalid argument 2: The Number of points must be greater than 10." << endl;
        exit(1);
    }

    vector<vector<double>> points;
    for (int i=0; i<dim; i++)
        points.push_back(Utils::createLejaSequence(size));
    Utils::storeLejaSequenceInFile(points);
    
    return 0;
}
