#!/bin/sh

make clean
make -j8

help()
{
  echo ""
  echo " To execute the script, you need to pass one argument:"
  echo "   - LEJA  : to see the leja sequence."
  echo "             Specify the dimension (2 or 3) and the number of points in one direction."
  echo "   - AI    : for AI algorithm."
  echo "   - PATH  : to see the path generated by AI algorithm."
  echo "   - PLOT  : to see the interpolation results and the progression of the algoritm."
  echo "   - ERROR : to see the interpolation error."
  echo ""
}

showAllArgsDetails()
{
  echo ""
  echo " 4 arguments are required:"
  echo "   - arg 1 : Space dimension [$2]"
  echo "   - arg 2 : Number of test points [$3]"
  echo "   - arg 3 : Method [$4]"
  echo "   - arg 4 : Number of iteration in AI algorithm [$5]"
  echo ""
}

show3ArgsDetails()
{
  echo ""
  echo " 4 arguments are required:"
  echo "   - arg 1 : Space dimension [$2]"
  echo "   - arg 2 : Number of test points [$3]"
  echo "   - arg 3 : Number of iteration in AI algorithm [$5]"
  echo ""
}

show2Details()
{
  echo ""
  echo " 2 arguments are required:"
  echo "   - arg 1 : Method [$2]"
  echo "   - arg 2 : Number of iteration in AI algorithm [$3]"
  echo ""
}

if [ "$1" = "LEJA" ]
then
    echo " - Leja sequence: "
    cd python
    python3.5 -W ignore plot_leja_sequence.py $2 $3

elif [ "$1" = "AI" ]
then
    showAllArgsDetails
    if [ $# != 5 ]
    then echo "Invalid number of arguments"
    else
        if [ $4 = 0 ]
        then ./bin/TestLagrangeInterpolation $2 $3 $5 0 0
      else ./bin/TestPiecewiseInterpolation $2 $3 $4 $5 0 0
        fi
    fi

elif [ "$1" = "MIXTE" ]
then
    show3ArgsDetails
    if [ $# != 4 ]
    then echo "Invalid number of arguments"
    else ./bin/TestMixedInterpolation $2 $3 $4 0 0
    fi

elif [ "$1" = "PATH" ]
then
    showAllArgsDetails
    if [ $# != 5 ]
    then echo "Invalid number of arguments"
    else
        if [ $4 = 0 ]
        then ./bin/TestLagrangeInterpolation $2 $3 $5 1 0
        else ./bin/TestPiecewiseInterpolation $2 $3 $4 $5 1 0
        fi
        cd python
        python3.5 -W ignore plot_path.py $2 $4
    fi

elif [ "$1" = "PLOT" ]
then
    show2Details
    if [ $# != 3 ]
    then echo "Invalid number of arguments"
    else
        if [ $2 = 0 ]
        then ./bin/TestLagrangeInterpolation 1 1000 $3 1 0
        else ./bin/TestPiecewiseInterpolation 1 1000 $2 $3 1 0
        fi
        cd python
        python3.5 -W ignore plot_interpolation_progression.py
    fi

elif [ "$1" = "ERROR" ]
then
    showAllArgsDetails
    if [ $# != 5 ]
    then echo "Invalid number of arguments"
    else
        if [ $4 = 0 ]
        then ./bin/TestLagrangeInterpolation $2 $3 $5 0 1
        else ./bin/TestPiecewiseInterpolation $2 $3 $4 $5 0 1
        fi
        cd python
        python3.5 -W ignore plot_error.py
    fi
else
    help
    exit
fi
