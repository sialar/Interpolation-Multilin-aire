#include <python3.5/Python.h>
#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/Tucker/TuckerApproximation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    string core = "MOX";
    vector<string>  crossSections;
    crossSections.push_back("macro_totale0");
    /*
    crossSections.push_back("macro_totale1");
    crossSections.push_back("macro_absorption0");
    crossSections.push_back("macro_absorption1");
    crossSections.push_back("macro_scattering000");
    crossSections.push_back("macro_scattering001");
    crossSections.push_back("macro_scattering010");
    crossSections.push_back("macro_scattering011");
    crossSections.push_back("macro_nu*fission0");
    crossSections.push_back("macro_nu*fission1");
    crossSections.push_back("macro_fission0");
    crossSections.push_back("macro_fission1");
    */
    TuckerApproximation TA(core,crossSections);
    MultiVariatePoint<double> point(5,0,0);
    point(0) = 40000;
    point(1) = 1000;
    point(2) = 0.5;
    point(3) = 0;
    point(4) = 0.5;
    double value_T;
    vector<double> values;
    for (string csName : TA.listOfCrossSectionNames)
    {
        value_T = TA.evaluate(point,csName);
        values.push_back(value_T);
        cout << value_T << endl;
    }

    return 0;
}
