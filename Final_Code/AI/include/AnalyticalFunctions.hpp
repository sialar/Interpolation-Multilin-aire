#ifndef ANALYTICALFUNCTIONS
#define ANALYTICALFUNCTIONS

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <stdio.h>

#include "MultiVariatePoint.hpp"
#include "Functions.hpp"

using namespace std;

class AnalyticalFunctions : public Functions
{
    private:
    	int m_d, m_n;

    public:
      	~AnalyticalFunctions() {};
      	AnalyticalFunctions() {};
      	AnalyticalFunctions(int d, int n);

      	vector<double> evaluate(MultiVariatePoint<double> x);
};

typedef std::shared_ptr<AnalyticalFunctions> AnalyticalFunctionsPtr;

#endif
