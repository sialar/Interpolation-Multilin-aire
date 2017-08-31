#!/bin/sh

help()
{
  echo ""
  echo " 1 argument est requis:"
  echo "   - FA : pour interpoler une fonction analytique"
  echo "   - DR : pour interpoler des données réelles"
  echo "   - doc : pour générer la documentation du code"
  echo ""
}

helpFA()
{
  echo ""
  echo " 4 arguments sont requis:"
  echo "   - arg 1  : Dimension de l'espace de départ de l'interpolé f"
  echo "   - arg 2  : Dimension de l'espace d'arrivé de l'interpolé f"
  echo "   - arg 3  : Nombre d'itérations dans l'algorithme AI"
  echo "   - arg 4  : Nombre de points de test pour l'évaluation de la méthode"
  echo ""
}

helpDR()
{
  echo ""
  echo " 3 arguments sont requis:"
  echo "   - arg 1  : Dimension de l'espace de départ de l'interpolé f (doit correspondre à la dimension spécifié dans le fichier de données \"AI/data/input.dat\")"
  echo "   - arg 2  : Dimension de l'espace d'arrivé de l'interpolé f (doit correspondre à la dimension spécifié dans le fichier de données \"AI/data/input.dat\")"
  echo "   - arg 3  : Nombre d'itérations dans l'algorithme AI"
  echo ""
}


if [ "$1" = "FA" ]
then
    helpFA
    shift
    AI/bin/TestWithAnalyticalFunction $*

elif [ "$1" = "DR" ]
then
    helpDR
    shift
    AI/bin/TestWithRealFunction $*

elif [ "$1" = "PATH" ]
then
    shift
    cd AI/python
    python3.5 -W ignore plot_path.py $1

elif [ "$1" = "PROGRESS" ]
then
    shift
    cd AI/python
    python3.5 -W ignore plot_interpolation_progression.py

elif [ "$1" = "POINTS" ]
then
    shift
    cd AI/python
    python3.5 -W ignore plot_interpolation_points_1D.py $1

elif [ "$1" = "DOC" ]
then
    cd AI/doc
    firefox html/index.html

else
    help
    exit
fi
