#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include <limits>
#include <stdio.h>
#include "../include/Utils.hpp"
#include "../include/Functions.hpp"
#include "../include/Tucker/LagrangePolynomial.hpp"
#include "../include/Tucker/TuckerApproximation.hpp"


using namespace std;



int main( int argc, char* argv[] )
{
    string s1 = "-1.45647840e-03+0.j";
    string s2 = "7.74139602e-03+0.j";

    cout << Utils::getRealPart(s1) << endl;
    cout << Utils::getRealPart(s2) << endl;

    return 0;
}
