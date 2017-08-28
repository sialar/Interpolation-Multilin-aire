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

  protected:
    	// m_d : dimension de l'éspace de départ de l'interpolé
    	// m_n : dimension de l'éspace d'arrivé de l'interpolé
    	int m_d, m_n;
      // Domaine de définition de chaque variable [-1,1] par défaut
      vector<vector<double>> m_parametersDomain;

  public:
    	~Functions() {};
    	Functions() {};
    	Functions(int d, int n);

      int getD() { return m_d;};
      int getN() { return m_n;};

	     // Implémente une fonction analytique
      virtual vector<double> evaluate(MultiVariatePoint<double> x) = 0;
      void displayParametersDomain();
};

typedef std::shared_ptr<Functions> FunctionsPtr;

#endif
