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
	// m_d : dimension de l'éspace de départ de l'interpolé
	// m_n : dimension de l'éspace d'arrivé de l'interpolé
	int m_d, m_n;


    public:
      	~Functions() {};
      	Functions() {};
      	Functions(int d, int n);

        int getD() { return m_d;};
        int getN() { return m_n;};

	       // Implémente une fonction analytique
        virtual vector<double> evaluate(MultiVariatePoint<double> x) = 0;
};

typedef std::shared_ptr<Functions> FunctionsPtr;

#endif
