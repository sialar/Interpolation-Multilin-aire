#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include <limits>
#include <stdio.h>
#include "../include/Tucker/TuckerApproximation.hpp"
#include "../include/Utils.hpp"

using namespace std;

void displayDomain(vector<vector<double>> domain)
{
    int m_d = 5;
    cout << " - Real domain of the cross section function : ";
    for (int i=0; i<m_d-1; i++)
      cout << "[" << domain[i][0] << "," << domain[i][1] << "] x ";
    cout << "[" << domain[m_d-1][0] << "," << domain[m_d-1][1] << "]" << endl;
}

int main( int argc, char* argv[] )
{
    MultiVariatePoint<double> x(5,0,0);
    x(0) = -1;
    x(1) = -1;
    x(2) = -1;
    x(3) = -1;
    x(4) = -1;

    vector<vector<double>> realDomain;
    vector<double> temp(2);
    temp[0] = 0;
    temp[1] = 80000;
    realDomain.push_back(temp);
    temp[0] = 20;
    temp[1] = 1800;
    realDomain.push_back(temp);
    temp[0] = 0.4;
    temp[1] = 1;
    realDomain.push_back(temp);
    temp[0] = 0;
    temp[1] = 1.74835e-05;
    realDomain.push_back(temp);
    temp[0] = 1e-06;
    temp[1] = 2;
    realDomain.push_back(temp);
    displayDomain(realDomain);

    cout << x << endl;
    for (int i=0; i<5; i++)
        x(i) = Utils::convertToFunctionDomain(realDomain[i][0],realDomain[i][1],x(i));
    cout << x << endl;
    for (int i=0; i<5; i++)
        x(i) = Utils::adaptCoordsToFunctionDomain(realDomain[i][0],realDomain[i][1],x(i));
    cout << x << endl;

    return 0;
}
