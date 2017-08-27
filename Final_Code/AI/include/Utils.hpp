#ifndef UTILS
#define UTILS

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>

#include "MultiVariatePoint.hpp"

using namespace std;

class Utils
{
    private:
	// Intervient dans la construction des points de Leja
        static vector<double> m_1dGrid;

    public:
   	
	// Fonctions d'affichage
        static void separateur();
        static void displayValues(vector<double> v);

	// Fonctions intervenants dans le calcul des erreurs
        static double min_elt(vector<double> v);
        static double max_elt(vector<double> v);
	static double norm(vector<double> x, int p);
        static double maxAbsValue(vector<double> vec);
        static vector<double> diff(vector<double> x,vector<double> y);
        static double computeMseError(vector<double> f,vector<double> f_tilde);
	static double relativeError(vector<double> f,vector<double> f_tilde);

	// Opérations de changement de variables
        static double adaptCoordsToFunctionDomain(double a, double b, double x);
        static double convertToFunctionDomain(double a, double b, double x);

	// Opérations sur des chaines de caractères
        static string replace(string strs, string str_old, string str_new);
        static string eraseExtraSpaces(string strs);
        static bool strInVector(string required, vector<string> vec);
        static string vector2str(vector<double> x);
        static vector<double> str2vector(string line);

	// Fonctions pour construire les séquences de points
        static double randomValue(double a, double b);
        static MultiVariatePoint<double> createRandomMultiVariatePoint(int d);
        static vector<double> createUniformSequence(int nbPoints);
        static vector<double> createChebychevSequence(int nbPoints);
        static vector<double> createLejaSequence(int nbPoints);
        static vector<double> loadLejaSequenceFromFile(int length);
        static void storeLejaSequenceInFile(int length);
        static void binaryDecomposition(int number, vector<double>& binary_decomp);
        static bool isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold);
        static double computeNewLejaPointFromSequence(vector<double> seq);
	
	// Comparaison de points multivariés
        static bool equals(MultiVariatePoint<string> nu1, MultiVariatePoint<string> nu2);
};

#endif
