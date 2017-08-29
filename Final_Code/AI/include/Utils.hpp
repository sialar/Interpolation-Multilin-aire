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

/**
 *  \file Utils.hpp
 *  \brief Classe statique qui implémente certaines opérations et fonctions utiles.
 *  \author SIALA Rafik
 *  \date 08/16
*/
class Utils
{
    private:
	     // Intervient dans la construction des points de Leja
       // La recherche du point suivant dans la construction des points de Leja se fait parmi
       // l'ensemble de points m_1dGridœ
       static vector<double> m_1dGrid;

    public:
        // Paramètre utilisé pour la precision d'affichage des doubles
        static double m_precision;

/******************************************************************************/
/************************ Fonctions d'affichage *******************************/
        /**
         *  Fonction pour afficher les élément d'un vecteur
         *  \param v : vecteur à afficher
        */
        static void displayValues(vector<double> v);

/******************************************************************************/
/****************** Opérations sur les vecteur de double **********************/
        /**
         *  Chercher le mix d'un vecteur
         *  \param v : vecteur de double
         *  \return element min de v
        */
        static double min_elt(vector<double> v);
        /**
         *  Chercher le max d'un vecteur
         *  \param v : vecteur de double
         *  \return element max de v
        */
        static double max_elt(vector<double> v);
        /**
         *  Calculer la norme d'un vecteur v
         *  si p =  0 --> norme infini
         *  si p != 0 --> norme 2
         *  \param v : vecteur de double
         *  \return norme de x
        */
        static double norm(vector<double> x, int p);
        /**
         *  Calculer la difference entre 2 vecteurs x et y
         *  \param x : vecteur de double
         *  \param y : vecteur de double
         *  \return x - y = (x_i - y_i)_i
        */
        static vector<double> diff(vector<double> x,vector<double> y);

/******************************************************************************/
/**************** Opérations de changement de variables ***********************/
        /**
         *  Changement de variable de [a,b] à [-1,1]
         *  \param a : borne inf du ségment d'entrée
         *  \param b : borne sup du ségment d'entrée
         *  \param x : point à transformer
         *  \return x dans [-1,1]
        */
        static double convertToDefaultDomain(double a, double b, double x);
        /**
         *  Changement de variable de [-1,1] à [a,b]
         *  \param a : borne inf du ségment de sortie
         *  \param b : borne sup du ségment de sortie
         *  \param x : point à transformer
         *  \return x dans [a,b]
        */
        static double convertToFunctionDomain(double a, double b, double x);

/******************************************************************************/
/**************** Opérations sur des chaines de caractères ********************/
        /**
         *  Remplacer dans strs les sous-chaines str_old par str_new
         *  \param strs : chaine de caractère
         *  \param str_old : sous-chaine remplacé
         *  \param new_old : nouvelle sous-chaine
         *  \return nouveau strs
        */
        static string replace(string strs, string str_old, string str_new);
        /**
         *  Supprimer les espaces de séparation inutiles dans strs
         *  En sortie strs est une chaine de caractères contenant des réels séparés par un seul espace
         *  \param strs : chaine de caractère
         *  \return strs sans espaces inutiles
        */
        static string eraseExtraSpaces(string strs);
        /**
         *  Vérifier si la chaine required existe dans le vecteur de chaine vec
         *  \param required : chaine de caractère recherchée
         *  \param vec : vecteur de chaines de caractère
         *  \return True si et seulement si requirede st un element de vec
        */
        static bool strInVector(string required, vector<string> vec);
        /**
         *  Transformer le vecteur x en chaine de caractère en séparant les double par un espace
         *  \param x : vecteur de double
         *  \return x sous forme de chaine de caractere
        */
        static string vector2str(vector<double> x);
        /**
         *  Transformer la chaine de caractère line en vecteur de double
         *  \param line : chaine de caractere à transformer
         *  \return vecteur de double
        */
        static vector<double> str2vector(string line);

/******************************************************************************/
/******** Foncions utiles pour la création des points d'interpolation *********/
        /**
         *  Retourner une valeur entre a et b aléatoirement
         *  \param a : borne inf
         *  \param b : borne sup
         *  \return valeur entre a et b
        */
        static double randomValue(double a, double b);
        /**
         *  Retourner un point multivarié dans [-1,1]^d d'une manière aléatoire
         *  \param d : taille du point
         *  \return point multivarié aléatoire
        */
        static MultiVariatePoint<double> createRandomMultiVariatePoint(int d);
        /**
         *  Créer une sequence uniforme de nbPoints points 1d dans [-1,1]
         *  \param nbPoints : taille de la séquence
         *  \return séquence uniforme de points
        */
        static vector<double> createUniformSequence(int nbPoints);
        /**
         *  Créer la sequence de Chebychev de nbPoints points 1d
         *  \param nbPoints : taille de la séquence
         *  \return séquence de Chebychev
        */
        static vector<double> createChebychevSequence(int nbPoints);
        /**
         *  Créer la sequence de Leja de nbPoints points 1d
         *  \param nbPoints : taille de la séquence
         *  \return séquence de Leja
        */
        static vector<double> createLejaSequence(int nbPoints);
        /**
         *  Charger la sequence de Leja de nbPoints points 1d depuis le fichier
         *  data/leja_sequence.dat
         *  \param nbPoints : taille de la séquence
         *  \return séquence de Leja
        */
        static vector<double> loadLejaSequenceFromFile(int nbPoints);
        /**
         *  Créer puis sauvegarder dans le fichier data/leja_sequence.dat la sequence de Leja
         *  de nbPoints points 1d
         *  \param nbPoints : taille de la séquence
        */
        static void storeLejaSequenceInFile(int nbPoints);
        /**
         *  Calculer l'expansion binaire d'un entier number
         *  \param number : Entier à décomposer
         *  \param binary_decomp : expansion binaire de number sous forme de vecteur
        */
        static void binaryExpansion(int number, vector<double>& binary_decomp);
        /**
         *  Verifier si le point y est trés proche d'un point de la séquence seq
         *  \param y : point concidéré
         *  \param seq : sequence de points 1d
         *  \param threshold : seuil de comparaison
         *  \return True si la distance de y à un point de seq est inférieure à threshold
        */
        static bool isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold);
        /**
         *  Calculer le nouveau point de Leja dans la liste de points seq
         *  \param seq : sequence courante de Leja
         *  \return point de Leja suivant dans la sequence seq
        */
        static double computeNewLejaPointFromSequence(vector<double> seq);

/******************************************************************************/
/************** Foncions de comparaison de points multivariés *****************/
        /**
         *  Test d'égalité entre les points multivariés nu1 et nu2
         *  \param nu1 : point multivarié
         *  \param nu2 : point multivarié
         *  \return True si et seulement si nu1 = nu2
        */
        static bool equals(MultiVariatePoint<string> nu1, MultiVariatePoint<string> nu2);

/******************************************************************************/
/********** Foncions utils pour la lecture des entrées du programme ***********/
        /**
        *  Choisir la dimension d de l'espace de départ de f
        *  \param argc : nombre d'arguments dans l'execution (arguments du main)
        *  \param argv : liste des arguments dans l'execution (arguments du main)
        *  \param argNum : indice de l'argument concidéré
        *  \return dimension d
        */
        static int chooseDimensionD(int argc, char* argv[], int argNum);
        /**
        *  Choisir la dimension n de l'espace d'arrivé de f
        *  \param argc : nombre d'arguments dans l'execution (arguments du main)
        *  \param argv : liste des arguments dans l'execution (arguments du main)
        *  \param argNum : indice de l'argument concidéré
        *  \return dimension n
        */
        static int chooseDimensionN(int argc, char* argv[], int argNum);
        /**
        *  Choisir le nombre de points pour tester l'interpolant
        *  \param argc : nombre d'arguments dans l'execution (arguments du main)
        *  \param argv : liste des arguments dans l'execution (arguments du main)
        *  \param argNum : indice de l'argument concidéré
        *  \return nombre de points de test
        */
        static int chooseNbTestPoints(int argc, char* argv[], int argNum);
        /**
        *  Choisir le nombre d'itérations (ou nombre de points d'interpolation) dans l'algorithme AI
        *  \param argc : nombre d'arguments dans l'execution (arguments du main)
        *  \param argv : liste des arguments dans l'execution (arguments du main)
        *  \param argNum : indice de l'argument concidéré
        *  \return nombre d'itérations
        */
        static int chooseMaxIteration(int argc, char* argv[], int argNum);
        /**
        *  Choisir la version de l'algorithme AI à utiliser dans chaque direction.
        *  Utile pour la version MixedInterpolation
        *  \param dim : dimension d de l'espace de depart de f
        *  \return multi-indice, chaque indice correspond à une version de AI
        */
        static MultiVariatePoint<int> chooseMethods(int dim);
};

#endif
