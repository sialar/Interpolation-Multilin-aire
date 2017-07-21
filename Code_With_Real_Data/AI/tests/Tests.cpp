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
    string s1 = "8.96460098e-01 -0.00000000e+00j";
    string s2 = "-2.35105019e-02-0.11393989j";

    cout << stod(Utils::getRealPart(s1) )<< endl;
    cout << stod(Utils::getRealPart(s2)) << endl;
    cout << stod(s2) << endl;

    return 0;
}
