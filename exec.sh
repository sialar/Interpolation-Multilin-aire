#!/bin/sh

#make clean
make -j8

help()
{
  echo ""
  echo " To launch the script, you need to choose one argument from the list below:"
  echo "   - LEJA   : to see the leja sequence."
  echo "   - AI     : to exectue the Adaptative Interpolation algorithm."
  echo "   - PATH   : to see the interpolation points generated by AI algorithm."
  echo "   - PLOT   : to see the progression of the approximation."
  echo "   - ERROR  : to see the interpolation error."
  echo "   - MIX    : to use different methods in different directions."
  echo "   - AUTO   : to let the algorithm choose the optimal method in each direction."
  echo "   - DIFF   : to see the interpolation results when using interpolation points obtained with different function."
  echo "   - TUCKER : to see results obtained with Tucker method."
  echo ""
}

showAllArgsDetails()
{
  echo ""
  echo " 5 arguments are required:"
  echo "   - arg 1 : Space dimension D [$2]"
  echo "   - arg 2 : Space dimension N [$3]"
  echo "   - arg 3 : Number of test points [$4]"
  echo "   - arg 4 : Method [$5]"
  echo "   - arg 5 : Number of iteration in AI algorithm [$6]"
  echo ""
}

showErrorArgsDetails()
{
  echo ""
  echo " 5 arguments are required:"
  echo "   - arg 1 : Space dimension D [$2]"
  echo "   - arg 2 : Space dimension N [$3]"
  echo "   - arg 3 : Number of test points [$4]"
  echo "   - arg 3 : Method [$5] (ALL : to see all methods results)"
  echo "   - arg 4 : Number of iteration in AI algorithm [$6]"
  echo ""
}

show4ArgsDetails()
{
  echo ""
  echo " 4 arguments are required:"
  echo "   - arg 1 : Space dimension D [$2]"
  echo "   - arg 2 : Space dimension N [$3]"
  echo "   - arg 3 : Number of test points [$4]"
  echo "   - arg 4 : Number of iteration in AI algorithm [$5]"
  echo ""
}

show2ArgsDetails()
{
  echo ""
  echo " 2 arguments are required:"
  echo "   - arg 1 : Method [$2]"
  echo "   - arg 2 : Number of iteration in AI algorithm [$3]"
  echo ""
}

show1ArgsDetails()
{
  echo ""
  echo " 1 argument is required:"
  echo "   - arg : Sample size of accuracy measurment [$2]"
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
    if [ $# != 6 ]
    then echo "Invalid number of arguments"
    else
        if [ $5 = 0 ]
        then ./bin/TestLagrangeInterpolation $2 $3 $4 $6 0 0
      else ./bin/TestPiecewiseInterpolation $2 $3 $4 $5 $6 0 0
        fi
    fi

elif [ "$1" = "MIX" ]
then
    show4ArgsDetails
    if [ $# != 5 ]
    then echo "Invalid number of arguments"
  else ./bin/TestMixedInterpolation $2 $3 $4 $5
    cd python
    python3.5 -W ignore plot_mixed_path.py $2
    fi

elif [ "$1" = "AUTO" ]
then
    show4ArgsDetails
    if [ $# != 5 ]
    then echo "Invalid number of arguments"
    else ./bin/TestAutoMixedInterpolation $2 $3 $4 $5
    cd python
    python3.5 -W ignore plot_mixed_path.py $2
    fi

elif [ "$1" = "PATH" ]
then
    showAllArgsDetails
    if [ $# != 6 ]
    then echo "Invalid number of arguments"
    else
        if [ $5 = 0 ]
        then ./bin/TestLagrangeInterpolation $2 $3 $4 $6 1 0
        else ./bin/TestPiecewiseInterpolation $2 $3 $4 $5 $6 1 0
        fi
        cd python
        python3.5 -W ignore plot_path.py $2 $5
    fi

elif [ "$1" = "PLOT" ]
then
    show2ArgsDetails
    if [ $# != 3 ]
    then echo "Invalid number of arguments"
    else
        if [ $2 = 0 ]
        then ./bin/TestLagrangeInterpolation 1 1 1000 $3 1 0
        else ./bin/TestPiecewiseInterpolation 1 1 1000 $2 $3 1 0
        fi
        cd python
        python3.5 -W ignore plot_interpolation_progression.py
    fi

elif [ "$1" = "DIFF" ]
then
    showAllArgsDetails
    if [ $# != 6 ]
    then echo "Invalid number of arguments"
  else ./bin/TestSameFunctionWithDifferentPaths $2 $3 $4 $5 $6
    cd python
    python3.5 -W ignore plot_2_paths.py $2 $5
    fi

elif [ "$1" = "ERROR" ]
then
    showErrorArgsDetails
    if [ $# != 6 ]
    then echo "Invalid number of arguments"
    else
        if [ $5 = 0 ]
        then ./bin/TestLagrangeInterpolation $2 $3 $4 $6 0 1
      elif [ "$5" = "ALL" ]
        then ./bin/TestErrors $2 $3 $4 $6
        else ./bin/TestPiecewiseInterpolation $2 $3 $4 $5 $6 0 1
        fi
        cd python
        python3.5 -W ignore plot_error.py
    fi

elif [ "$1" = "TUCKER" ]
then
    show1ArgsDetails
    if [ $# != 2 ]
    then echo "Invalid number of arguments"
    else cd tucker
         python testTuckerDecomposition_withGreedy.py $2
    fi
else
    help
    exit
fi
