#!/bin/bash

#make -j8

if [ $# == 0 ]
then
    echo ""
    echo " To execute the script, you need to pass one argument:"
    echo "   - LEJA : for leja sequence computation."
    echo "     Specify the dimension (2 or 3) and the number of points in one direction."
    echo "   - AI   : for AI algorithm."
    echo "     Interactive execution."
    echo ""
fi

if [ "$1" == "LEJA" ]
then
    echo " - Sequence de Leja"
    ./bin/TestLejaSequence $2 $3
    cd python
    python3.5 leja_sequence.py

elif [ "$1" == "AI" ]
then
    echo ""
    echo " - Interpolation in dimension (d>0) + AI algorithm:"
    ./bin/TestAlgoAI $2
    if [ "$2" == "-p" ]
    then
        cd python
        python3.5 progressive_plot.py
    fi
fi
