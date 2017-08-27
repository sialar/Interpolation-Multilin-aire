#ifndef REALDATAFUNCTIONS
#define REALDATAFUNCTIONS

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
#include "Utils.hpp"

using namespace std;

class RealDataFunctions : public Functions
{
    private:
      	int m_d, m_n;
        string m_filePath;


    public:
      	~RealDataFunctions() {};
      	RealDataFunctions() {};
      	RealDataFunctions(int d, int n, string file);

      	vector<double> evaluate(MultiVariatePoint<double> x);
};

typedef std::shared_ptr<RealDataFunctions> RealDataFunctionsPtr;

#endif
