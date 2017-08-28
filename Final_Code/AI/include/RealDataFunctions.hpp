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

#include "Tucker/TuckerApproximation.hpp"
#include "MultiVariatePoint.hpp"
#include "Functions.hpp"
#include "Utils.hpp"

using namespace std;

class RealDataFunctions : public Functions
{
    private:
      	int m_d, m_n;
        string m_filePath;
        vector<MultiVariatePoint<double>> m_refPoints;
        vector<vector<double>> m_exactValues;

        string m_coreName = "MOX"; // MOX ou UOX ou UOX-Gd
        string m_csName = "macro_totale0"; // macro_totale1 ou maco_fission1 ou ...
        TuckerApproximationPtr m_tuckerApprox;

    public:
      	~RealDataFunctions() {};
      	RealDataFunctions() {};
      	RealDataFunctions(int d, int n, string file);

      	vector<double> evaluate(MultiVariatePoint<double> x);
        void displayParametersDomain();
        void readDataFromFile();
};

typedef std::shared_ptr<RealDataFunctions> RealDataFunctionsPtr;

#endif
