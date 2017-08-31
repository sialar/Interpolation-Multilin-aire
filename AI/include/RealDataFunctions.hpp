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

/**
 *  \file RealDataFunctions.hpp
 *  \brief Classe dérivée de Functions: pour approcher des fonctions présentées sous forme de données réelles.
 *  \author SIALA Rafik
 *  \date 08/16
*/
class RealDataFunctions : public Functions
{
    /*
     * Il faut avoir:
     *  - Un fichier (eg. cross_section.dat) qui contient un ensemble de points et la valeur de f en chaque point.
     *    Les données doivent être organisées comme dans le fichier cross_section.dat
     *    Ce fichier constitue, en pratique, la seule information qu'on a sur f
     *  - Un code (doit être implémenté dans evaluate(x)) qui permet de donner la valeur "exacte" (de référence) de f
     *    au point x. x peut être différent de tous les points de référence dans le fichier m_fileName.
     *    En effet, ce qui est souvent le cas, quand l'algo AI construit d'une manière adaptative un point
     *    d'interpolation et utilise la valeur exacte en ce point.
     *    En pratique, ce code peut donner une valeur approximative (qu'on suppose exacte) de f en un point.
    */
    private:
        // m_d : dimension de l'éspace de départ de f
        // m_n : dimension de l'éspace d'arrivé de f
      	int m_d, m_n;
        // Nom du fichier contenant les données.
        // Les données comprennent les points de référence et la valeur de la fonction f
        // qu'on veut approcher en ces points. eg. cross_section.dat
        string m_fileName;

        // Ensemble des points de référence (lus à partir du fichier m_fileName)
        vector<MultiVariatePoint<double>> m_referencePts;
        // Ensemble des valeurs de f (lus à partir du fichier m_fileName)
        vector<vector<double>> m_exactValues;

        // Exemple d'utilisation
        // - m_fileName contient les valeur de la section efficace totale de type MOX.
        // - On essaie d'approcher la valeur obtenu par APOLLO2 (supposée exacte).
        // - On n'a pas le code APOLLO2. On va donc approcher la section efficace en utilisant un code
        // ( Modèle de décomposition de Tucker)
        // Cet attribut donne accés à la classe qui permet d'approcher les sections efficaces par la
        // méthode de décomposition de Tucker.
        TuckerApproximationPtr m_tuckerApprox;

    public:
      	~RealDataFunctions() {};
      	RealDataFunctions() {};
        /**
         *  Constructeur
         *  \param d: dimension de l'éspace de départ de f
         *  \param n: dimension de l'éspace de départ de f
         *  \param file: nom du fichier de données
        */
      	RealDataFunctions(int d, int n, string file);

        /**
         *  Retourner la liste des points de référence lus à partir du fichier m_fileName
         *  \return Ensemble des points de références
        */
        vector<MultiVariatePoint<double>> referencePts() { return m_referencePts; };
        /**
         *  Retourner la liste des valeurs de f aux points de référence
         *  \return Ensemble des valeurs de f (peuvent être des valeurs vectorielles)
        */
        vector<vector<double>> exactValues() { return m_exactValues; };

        /**
         *  Retourner la valeur "exacte" de f au point multivariée x
         *  \param x : points multivarié
         *  \return f(x) : valeur de f au point x
        */
        vector<double> evaluate(MultiVariatePoint<double> x);

        /**
         *  Lire les données à partir du fichier m_fileName et les stocker dans les attributs
         *  m_referencePts et m_exactValues. Ils serviront à construire puis à évaluer l'opérateur d'interpolation.
        */
        void readDataFromFile();
};

typedef std::shared_ptr<RealDataFunctions> RealDataFunctionsPtr;

#endif
