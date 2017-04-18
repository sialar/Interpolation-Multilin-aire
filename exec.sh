#!/bin/bash

if [ $# > 0 ]
then
    make clean
    make
    if [ $1 = "1d" ]
    then
        echo ""
        echo "Interpolation 1d:"
        echo ""
        ./bin/test1D $2 $3

    elif [ $1 = "2d" ]
    then
        echo ""
        echo "Interpolation 2d:"
        echo ""
        ./bin/test2D $2 $3 $4 $5 $6
        cd python
        python3.5 show_result.py

    elif [ $1 = "leja" ]
    then
        ./bin/testLejaSequence $2 $3
        echo "Affichage de $2 x $3 points de Leja"
        cd python
        python3.5 leja_sequence.py
    fi
fi
