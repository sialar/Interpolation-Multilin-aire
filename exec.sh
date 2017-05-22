#!/bin/bash

make clean
make -j8

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
    python3.5 -W ignore leja_sequence.py

elif [ "$1" == "AI" ]
then
    echo ""
    echo " 4 arguments are required:"
    echo "   - arg 1 : Space dimension [$2]"
    echo "   - arg 2 : Number of test points [$3]"
    echo "   - arg 3 : Method [$4]"
    echo "   - arg 4 : Number of iteration in AI algorithm [$5]"
    echo ""
    ./bin/TestAlgoAI $2 $3 $4 $5
elif [ "$1" == "PATH" ]
then
    echo ""
    ./bin/TestAlgoAI 2 1 0
    cd python
    python3.5 -W ignore progressive_plot.py
elif [ "$1" == "PLOT" ]
then
    echo ""
    ./bin/TestAlgoAI 1 10000 $2
    cd python
    python3.5 -W ignore plot_function.py $2
fi
