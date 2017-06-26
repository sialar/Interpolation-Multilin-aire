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
enum CoreType { MOX, UOX, UOX_Gd  };
enum ReactionType { ABS_0, ABS_1, FIS_0, FIS_1, NU_FIS_0, NU_FIS_1, \
                    SCAT_00, SCAT_01, SCAT_10, SCAT_11, TOT_0, TOT_1 };

class Functions
{
    private:
      int m_n;
      CoreType m_coreType;
      vector<ReactionType> m_reactionTypes;
      map<MultiVariatePoint<double>, vector<double>> m_values;

    public:
      ~Functions() {};
      Functions() {};
      Functions(int n, CoreType c, vector<ReactionType> vr);

      CoreType getCoreType() { return m_coreType; };
      void setCoreType(CoreType c) { m_coreType = c; };

      vector<ReactionType> getReactionType() { return m_reactionTypes; };
      void setReactionTypes(vector<ReactionType> vr);

      string getFunctionFileName(CoreType, ReactionType);
      vector<double> evaluate(MultiVariatePoint<double> x);
      void setFunctionsMap();
      void updateValues();
};

#endif
