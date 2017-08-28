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

/**
 *  \file AnalyticalFunctions.hpp
 *  \brief Classe dérivée de Functions: implémente des exemples de fonctions analytiques. Permet de valider la méthode.
 *  \author SIALA Rafik
 *  \date 08/16
*/
class AnalyticalFunctions : public Functions
{
    private:
      // m_d : dimension de l'éspace de départ de f
      // m_n : dimension de l'éspace d'arrivé de f
    	int m_d, m_n;

    public:
      	~AnalyticalFunctions() {};
      	AnalyticalFunctions() {};
        /**
         *  Constructeur
         *  \param d: dimension de l'éspace de départ de f
         *  \param n: dimension de l'éspace de départ de f
        */
      	AnalyticalFunctions(int d, int n);

        /**
         *  Retourner la valeur "exacte" de f au point multivariée x
         *  \param x : points multivarié
         *  \return f(x) : valeur de f au point x
        */
      	vector<double> evaluate(MultiVariatePoint<double> x);
};

typedef std::shared_ptr<AnalyticalFunctions> AnalyticalFunctionsPtr;

#endif
