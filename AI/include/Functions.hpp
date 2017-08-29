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

/**
 *  \file Functions.hpp
 *  \brief Classe de base abstraite modélisant les diffèrents types de fonctions qu'on veut approcher (fonctions analytiques ou données réelles)
 *  \author SIALA Rafik
 *  \date 08/16
*/
class Functions
{

  protected:
      // Fonction à approcher :  f: (x_1, .., x_{m_d}) --> (y_1, .., y_{m_n})
    	// m_d : dimension de l'éspace de départ de f
    	// m_n : dimension de l'éspace d'arrivé de f
    	int m_d, m_n;
      // Domaine de définition de chaque variable: [-1,1] par défaut
      vector<vector<double>> m_parametersDomain;

  public:
    	~Functions() {};
    	Functions() {};
      /**
        *  Constructeur
        *  \param d : dimension de l'éspace de départ de f
        *  \param n : dimension de l'éspace d'arrivé de f
      */
    	Functions(int d, int n);

      /**
        *  Accesseur à l'attribut m_d
        *  \return valeur de m_d
      */
      int getD() { return m_d;};
      /**
        *  Accesseur à l'attribut m_n
        *  \return valeur de m_n
      */
      int getN() { return m_n;};
      /**
        *  Accesseur à l'attribut m_parametersDomain
        *  \return domaine de définition de f
      */
      vector<vector<double>> parametersDomain() { return m_parametersDomain; };

      /**
        *  Méthode donner la valeur "exacte" de f en x
        *  \param x : point de calcul
        *  \return valeur de f en x
      */
      virtual vector<double> evaluate(MultiVariatePoint<double> x) = 0;
      /**
        *  Afficher le domaine de définition de f
      */
      void displayParametersDomain();
};

typedef std::shared_ptr<Functions> FunctionsPtr;

#endif
