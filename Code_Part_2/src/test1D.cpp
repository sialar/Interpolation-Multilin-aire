#include <iostream>
#include <string>
#include <time.h>
#include <algorithm>
#include "../include/Interpolation.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    if (argc != 2) return 0;

    int size = stoi(argv[1]);
    Interpolation* interp = new Interpolation(size);


    double testPoint = rand()/(double)RAND_MAX;
    double real = Interpolation::g(testPoint);
    double approx = interp->interpolation_iterative(testPoint,size,false);
    interp->displayPath();
    
    cout << endl << "Test point:" << testPoint << endl;
    cout << "Real value: g(" << testPoint << ") = " << real << endl;
    cout << "Aproximation value: ĝ(" << testPoint << ") = " << approx << endl;

    return 0;
}
