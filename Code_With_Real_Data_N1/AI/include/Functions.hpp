#ifndef FUNCTIONS
#define FUNCTIONS

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <stdio.h>

#include "MultiVariatePoint.hpp"
#include "Utils.hpp"

#define BUFSIZE 1024

using namespace std;

class Functions
{
    private:
      string m_directory;
      string m_coreType;
      string m_crossSection;

    public:
      ~Functions() {};
      Functions() {};
      Functions(string c, string cs);

      static vector<string> allCoreTypes;
      static vector<string> allCrossSectionType;

      string realDataDirPath() { return m_directory; };

      string getCoreType() { return m_coreType; };
      void setCoreType(string c);

      string getCrossSection() { return m_crossSection; };
      void setCrossSectionType(string cs);

      double evaluate(MultiVariatePoint<double> x);

      static void createFunctionsDataBase();
      static bool validCoreType(string c);
      static bool validCrossSection(string r);
};

typedef std::shared_ptr<Functions> FunctionsPtr;

#endif
