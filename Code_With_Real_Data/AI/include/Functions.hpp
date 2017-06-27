#ifndef FUNCTIONS
#define FUNCTIONS

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>

#include "MultiVariatePoint.hpp"
#include "Utils.hpp"

using namespace std;

class Functions
{
    private:
      int m_n;
      string m_directory;
      string m_coreType;
      vector<string> m_reactionTypes;


    public:
      ~Functions() {};
      Functions() {};
      Functions(int n, string c, vector<string> vr);

      static vector<string> allCoreTypes;
      static vector<string> allReactionTypes;

      string getCoreType() { return m_coreType; };
      void setCoreType(string c);

      vector<string> getReactionType() { return m_reactionTypes; };
      void setReactionTypes(vector<string> vr);
      void setAllReactionTypes();

      vector<double> evaluate(MultiVariatePoint<double> x);

      static void createFunctionsDataBase();
      static bool validCoreType(string c);
      static bool validReactionType(string r);
};

typedef std::shared_ptr<Functions> FunctionsPtr;

#endif
