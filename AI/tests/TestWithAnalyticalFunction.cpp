#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/MixedInterpolation.hpp"
#include "../include/LagrangeInterpolation.hpp"
#include "../include/AnalyticalFunctions.hpp"
#include "../include/Utils.hpp"

using namespace std;

// 3 arguments sont requis:
//   - arg 1  : Dimension de l'espace de départ de l'interpolé f (doit correspondre à la dimension spécifié dans le fichier de données \"AI/data/input.dat\")
//   - arg 2  : Dimension de l'espace d'arrivé de l'interpolé f (doit correspondre à la dimension spécifié dans le fichier de données \"AI/data/input.dat\")
//   - arg 3  : Nombre d'itérations dans l'algorithme AI
//   - arg 4  : Nombre de points de test pour l'évaluation de la méthode


int main( int argc, char* argv[] )
{
    // f : x = (x1, .., xd) --> (y1, .., yn)
    // Choix de la dimension d (manière interactive ou passage en argument)
    int d = Utils::chooseDimensionD(argc, argv, 1);
    // Choix de la dimension n (manière interactive ou passage en argument)
    int n = Utils::chooseDimensionN(argc, argv, 2);
    // Choix du nombre d'itérations (manière interactive ou passage en argument)
  	int nbIter = Utils::chooseMaxIteration(argc, argv, 3);
    // Choix du nombre de points de référence (de test) (manière interactive ou passage en argument)
  	int nbTestPts = Utils::chooseNbTestPoints(argc, argv, 4);
    // Choix de la version de l'algo AI à utiliser sur chaque direction
    // eg. d=2, methods = (i,j) --> utiliser la version i sur la direction 0 et j sur la direction 1
    // 3 versions possibles de l'algorithme d'interpolation adaptative:
    //  + 0 : Polynômes de Lagrange définis globalement + Points de Leja
    //  + 1 : Fonctions affines par morceaux + Construction des points par dichotomie
    //  + 2 : Fonctions quadratiques par morceaux + Construction des points par dichotomie
    MultiVariatePoint<int> methods = Utils::chooseMethods(d);
    // Construction de la fonction analytique à approcher
    AnalyticalFunctionsPtr f = make_shared<AnalyticalFunctions>(d,n);

    // Construction de l'opérateur d'interpolation en utilisant une combinaison des
    // différentes version (car utilisation de MixedInterpolation)
    MixedInterpolationPtr interp(new MixedInterpolation(f,nbIter,methods));
    // Construction d'une grille aléatoire de points de test
	  interp->setRandomTestPoints(nbTestPts);

    // Lancer l'algorithme d'interpolation adaptative
    interp->launchAIAlgo(false);
    // Caluler et stocker (dans output.dat) les résultats d'approximation sur la même grille construite aléatoirement
    interp->computeAIApproximationResults();
    interp->saveAll();

    // Afficher les résultats (precision, temps de calcul, nombre de points de calcul, ...) sur le terminal
    interp->displayAll();

    return 0;
}
