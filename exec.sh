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
    cd python
    python3.5 -W ignore plot_leja_sequence.py $2 $3

elif [ "$1" == "AI" ]
then
    echo ""
    echo " 4 arguments are required:"
    echo "   - arg 1 : Space dimension [$2]"
    echo "   - arg 2 : Number of test points [$3]"
    echo "   - arg 3 : Method [$4]"
    echo "   - arg 4 : Number of iteration in AI algorithm [$5]"
    echo ""
    if [ "$4" == "0" ]
    then
        ./bin/TestLagrangeInterpolation $2 $3 $5
    else
        ./bin/TestPiecewiseInterpolation $2 $3 $4 $5
    fi
elif [ "$1" == "PATH" ]
then
    echo ""
    ./bin/TestLagrangeInterpolation 2 1
    cd python
    python3.5 -W ignore plot_path.py

elif [ "$1" == "PLOT" ]
then
    if [ $# -le 3 ]
    then
        echo "Invalid number of arguments. Choose the method and the number of iteration"
    else
        if [ "$4" == "0" ]
        then
            ./bin/TestLagrangeInterpolation 1 1000 $5
        else
            ./bin/TestPiecewiseInterpolation 1 1000 $4 $5
        fi
    cd python
    python3.5 -W ignore plot_basis_functions.py
    fi
else
    echo ""
    echo " To execute the script, you need to pass one argument:"
    echo "   - LEJA : for leja sequence computation."
    echo "     Specify the dimension (2 or 3) and the number of points in one direction."
    echo "   - AI   : for AI algorithm."
    echo "     Interactive execution."
    echo ""
  fi
