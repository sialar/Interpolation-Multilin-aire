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

#include "Tucker/TuckerApproximation.hpp"
#include "MultiVariatePoint.hpp"
#include "Utils.hpp"

#define BUFSIZE 1024

using namespace std;

class Functions
{
    private:
      int m_n;
      string m_directory;
      string m_coreType;
      vector<string> m_crossSections;
      TuckerApproximationPtr m_tuckerProgram;


    public:
      ~Functions() {};
      Functions() {};
      Functions(string c, vector<string> cs);

      static vector<string> allCoreTypes;
      static vector<string> allCrossSectionType;

      TuckerApproximationPtr tuckerProgram() { return m_tuckerProgram; };
      void setTuckerProgram();

      string getCoreType() { return m_coreType; };
      void setCoreType(string c);

      vector<string> getReactionType() { return m_crossSections; };
      void setCrossSectionType(vector<string> cs);
      void setAllCrossSectionType();

      vector<double> evaluate(MultiVariatePoint<double> x);

      static void createFunctionsDataBase();
      static bool validCoreType(string c);
      static bool validReactionType(string r);
};

typedef std::shared_ptr<Functions> FunctionsPtr;

#endif