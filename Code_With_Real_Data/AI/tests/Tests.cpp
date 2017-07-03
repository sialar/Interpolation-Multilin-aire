#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include <limits>
#include <stdio.h>
#include "../include/Tucker/TuckerApproximation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    string s = "8.39912e-323 1	0.0	20.0	0.40000000596	0.0	9.99999997475e-07	0.0125141628087	0.0125141618773	0.0125141618773	-0.00502291352693	-0.00502291355499";
    cout << s << endl;
    string ss = Utils::eraseExtraSpaces(s);
    cout << ss << endl;
    return 0;
}
