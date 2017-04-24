#!/bin/bash

make clean
make -j8

if [ "$1" == "1d" ]
then
    echo ""
    echo "  ---> INTERPOLATION 1D:"
    ./bin/test1D $2 $3
    if [ "$4" == "-p" ]
    then
        cd python
        python3.5 show_result_1D.py
    fi
elif [ "$1" == "2d" ]
then
    echo ""
    echo "  ---> INTERPOLATION 2D:"
    ./bin/test2D $2 $3 $4 $5 $6
    if [ "$7" == "-p" ]
    then
        cd python
        python3.5 show_result_2D.py
    fi
elif [ "$1" == "leja" ]
then
    if [ $# == 3 ]
    then
        echo "   - Sequence de Leja 2D"
        ./bin/testLejaSequence $2 $3
    elif [ $# == 4 ]
    then
        echo "   - Sequence de Leja 3D"
        ./bin/testLejaSequence $2 $3 $4
    fi
    cd python
    python3.5 leja_sequence.py
elif [ "$1" == "AI" ]
then
    echo ""
    echo "  ---> INTERPOLATION 2D + ALGO AI:"
    ./bin/testAI $2 $3 $4 $5 $6
    if [ "$7" == "-p" ]
    then
      echo ""
      echo "  ---> AFFICHAGE DE LA PROGRESSION DU CHEMIN D'INDICES:"
      cd python
      python3.5 progressive_plot.py
    fi
elif [ "$1" == "ND_AI" ]
then
    echo ""
    echo "  ---> INTERPOLATION ND + ALGO AI:"
    ./bin/testND_AI $2 $3 $4 $5 $6 $7 $8 $9
fi
