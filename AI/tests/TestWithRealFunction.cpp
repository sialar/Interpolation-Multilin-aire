#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/MixedInterpolation.hpp"
#include "../include/LagrangeInterpolation.hpp"
#include "../include/RealDataFunctions.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    // f : x = (x1, .., xd) --> (y1, .., yn)
    // Choix de la dimension d (manière interactive ou passage en argument)
    int d = Utils::chooseDimensionD(argc, argv, 1);
    // Choix de la dimension n (manière interactive ou passage en argument)
  	int n = Utils::chooseDimensionN(argc, argv, 2);
    // Choix du nombre d'itérations (manière interactive ou passage en argument)
    int nbIter = Utils::chooseMaxIteration(argc, argv, 3);
    // Choix de la version de l'algo AI à utiliser sur chaque direction
    // eg. d=2, methods = (i,j) --> utiliser la version i sur la direction 0 et j sur la direction 1
    // 3 versions possibles de l'algorithme d'interpolation adaptative:
    //  + 0 : Polynômes de Lagrange définis globalement + Points de Leja
    //  + 1 : Fonctions affines par morceaux + Construction des points par dichotomie
    //  + 2 : Fonctions quadratiques par morceaux + Construction des points par dichotomie
    MultiVariatePoint<int> methods = Utils::chooseMethods(d);
    // Choix du fichier contenant les données réelles
    string fileName = "input.dat";
    // Construction de la fonction à approcher (données réelles stockées dans le fichier fileName)
    RealDataFunctionsPtr f = make_shared<RealDataFunctions>(d,n,fileName);

    // Construction de l'opérateur d'interpolation en utilisant une combinaison des
    // différentes version (car utilisation de MixedInterpolation)
    MixedInterpolationPtr interp(new MixedInterpolation(f,nbIter,methods));
    // Construction de la grille des points de test en utilisant les points lus dans le fichier fileName
	  interp->setTestPoints(f->referencePts());

    // Lancer l'algorithme d'interpolation adaptative
    interp->launchAIAlgo(false);
    // Caluler et stocker (dans output.dat) les résultats d'approximation sur la même grille
    // de référence présene dans le fichier fileName
    interp->computeAIApproximationResults();
    interp->saveAIApproximationResults();

    // Afficher les résultats (precision, temps de calcul, nombre de points de calcul, ...) sur le terminal
    interp->displayAll();

    return 0;
}
