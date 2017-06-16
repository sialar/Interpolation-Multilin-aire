#!/bin/sh


PROJECT_PATH="/home/sialar/Stage/LaboJ_LLions/Code"
cd $PROJECT_PATH

#make clean
make -j8

help()
{
  echo ""
  echo " To launch the script, you need to choose one argument from the list below:"
  echo "   - LEJA     : to see the leja sequence."
  echo "   - AI       : to exectue the Adaptative Interpolation algorithm."
  echo "   - PATH     : to see the interpolation points generated by AI algorithm."
  echo "   - PLOT     : to see the progression of the approximation."
  echo "   - ERROR    : to see the interpolation error."
  echo "   - AUTO     : to let the algorithm choose the optimal method in each direction."
  echo "   - DIFF     : to see the interpolation results when using interpolation points obtained with different function."
  echo "   - TUCKER   : to see results obtained with Tucker method."
  echo "   - COMP     : to compare results of both methods (Tucker and AI algo)."
  echo ""
}

showAIArgsDetails()
{
  echo ""
  echo " 6 arguments are required:"
  echo "   - arg 1 : Space dimension D [$2]"
  echo "   - arg 2 : Space dimension N [$3]"
  echo "   - arg 3 : Number of test points [$4]"
  echo "   - arg 4 : Method (0, 1, 2 or MIX) [$5]"
  echo "   - arg 5 : Number of iteration in AI algorithm [$6]"
  echo "   - arg 6 : Function to interpolate [$7]"
  echo ""
}

showAUTOArgsDetails()
{
  echo ""
  echo " 5 arguments are required:"
  echo "   - arg 1 : Space dimension D [$2]"
  echo "   - arg 2 : Space dimension N [$3]"
  echo "   - arg 3 : Number of test points [$4]"
  echo "   - arg 4 : Number of iteration in AI algorithm [$5]"
  echo "   - arg 6 : Function to interpolate [$6]"
  echo ""
}

showERRORArgsDetails()
{
  echo ""
  echo " 6 arguments are required:"
  echo "   - arg 1 : Space dimension D [$2]"
  echo "   - arg 2 : Space dimension N [$3]"
  echo "   - arg 3 : Number of test points [$4]"
  echo "   - arg 4 : Number of iteration in AI algorithm [$6]"
  echo "   - arg 5 : Function to interpolate [$7]"
  echo "   - arg 6 : The starting iteration in plot [$6]"
  echo ""
}

showDIFFArgsDetails()
{
  echo ""
  echo " 6 arguments are required:"
  echo "   - arg 1 : Space dimension D [$2]"
  echo "   - arg 2 : Space dimension N [$3]"
  echo "   - arg 3 : Number of test points [$4]"
  echo "   - arg 4 : Method [$5]"
  echo "   - arg 5 : Number of iteration in AI algorithm [$6]"
  echo "   - arg 6 : Function to interpolate [$7]"
  echo ""
}

showCOMPArgsDetails()
{
  echo ""
  echo " 4 arguments are required:"
  echo "   - arg 1 : Mumber of test points [$2]"
  echo "   - arg 2 : Method (0, 1, 2 or MIX) [$3]"
  echo "   - arg 3 : Number of iteration in AI algorithm [$4]"
  echo "   - arg 4 : Function to interpolate [$5]"
  echo ""
}

showPLOTArgsDetails()
{
  echo ""
  echo " 3 arguments are required:"
  echo "   - arg 1 : Method [$2]"
  echo "   - arg 2 : Number of iteration in AI algorithm [$3]"
  echo "   - arg 3 : Function to interpolate [$4]"
  echo ""
}

showTUCKERArgsDetails()
{
  echo ""
  echo " 2 argument is required:"
  echo "   - arg 1 : Number of test points [$2]"
  echo "   - arg 2 : Function to interpolate [$3]"
  echo ""
}

if [ "$1" = "LEJA" ]
then
    echo " - Leja sequence: "
    cd python
    python3.5 -W ignore plot_leja_sequence.py $2 $3

elif [ "$1" = "AI" ]
then
    showAIArgsDetails
    if [ $# != 7 ]
    then echo "Invalid number of arguments"
    else
        if [ $5 = 0 ]
        then "$PROJECT_PATH"/bin/TestLagrangeInterpolation $2 $3 $4 $6 $7 0
        elif [ "$5" = "MIX" ]
        then "$PROJECT_PATH"/bin/TestMixedInterpolation $2 $3 $4 $6 $7 0
        else "$PROJECT_PATH"/bin/TestPiecewiseInterpolation $2 $3 $4 $5 $6 $7 0
        fi
    fi

elif [ "$1" = "AUTO" ]
then
    showAUTOArgsDetails
    if [ $# != 6 ]
    then echo "Invalid number of arguments"
    else "$PROJECT_PATH"/bin/TestAutoMixedInterpolation $2 $3 $4 $5 $6 0
    fi

elif [ "$1" = "PATH" ]
then
    showAIArgsDetails
    if [ $# != 7 ]
    then echo "Invalid number of arguments"
    elif [ $2 != 2 ] && [ $2 != 3 ]
    then echo "Space dimension (arg 1) must be 2 or 3"
    else
        if [ $5 = 0 ]
        then "$PROJECT_PATH"/bin/TestLagrangeInterpolation $2 $3 $4 $6 $7 1
        elif [ "$5" = "MIX" ]
        then "$PROJECT_PATH"/bin/TestMixedInterpolation $2 $3 $4 $6 $7 1
        else "$PROJECT_PATH"/bin/TestPiecewiseInterpolation $2 $3 $4 $5 $6 $7 1
        fi
        cd python
        if [ "$5" = "MIX" ]
        then python3.5 -W ignore plot_mixed_path.py $2
        else python3.5 -W ignore plot_path.py $2 $5
        fi
    fi

elif [ "$1" = "PLOT" ]
then
    showPLOTArgsDetails
    if [ $# != 4 ]
    then echo "Invalid number of arguments"
    else
        if [ $2 = 0 ]
        then "$PROJECT_PATH"/bin/TestLagrangeInterpolation 1 1 1000 $3 $4 1
        else "$PROJECT_PATH"/bin/TestPiecewiseInterpolation 1 1 1000 $2 $3 $4 1
        fi
        cd python
        python3.5 -W ignore plot_interpolation_progression.py
    fi

elif [ "$1" = "DIFF" ]
then
    showDIFFArgsDetails
    if [ $# != 7 ]
    then echo "Invalid number of arguments"
    else "$PROJECT_PATH"/bin/TestSameFunctionWithDifferentPaths $2 $3 $4 $5 $6 $7
    cd python
    python3.5 -W ignore plot_2_paths.py $2 $5
    fi

elif [ "$1" = "ERROR" ]
then
    showERRORArgsDetails
    if [ $# != 7 ]
    then echo "Invalid number of arguments"
    else "$PROJECT_PATH"/bin/TestErrors $2 $3 $4 $5 $6 $7
    cd python
    python3.5 -W ignore plot_error.py $6
  fi

elif [ "$1" = "TUCKER" ]
then
    showTUCKERArgsDetails
    if [ $# != 3 ]
    then echo "Invalid number of arguments"
    else cd tucker
         python testTuckerDecomposition_withGreedy.py $2 $3
    fi

elif [ "$1" = "COMP" ]
then
    showCOMPArgsDetails
    if [ $# != 5 ]
    then echo "Invalid number of arguments"
    else
        echo "\t\t\t\t\t\t\t\t\t-------------------------------------------"
        echo "\t\t\t\t\t\t\t\t\t-- Using Adaptative Interpolation method --"
        echo "\t\t\t\t\t\t\t\t\t-------------------------------------------"
        if [ $3 = 0 ]
        then "$PROJECT_PATH"/bin/TestLagrangeInterpolation 3 1 $2 $4 $5 0
        elif [ "$3" = "MIX" ]
        then "$PROJECT_PATH"/bin/TestMixedInterpolation 3 1 $2 $4 $5 0
        else "$PROJECT_PATH"/bin/TestPiecewiseInterpolation 3 1 $2 $3 $4 $5 0
        fi
        echo "\n"
        echo "\t\t\t\t\t\t\t\t\t\t-------------------------"
        echo "\t\t\t\t\t\t\t\t\t\t-- Using Tucker method --"
        echo "\t\t\t\t\t\t\t\t\t\t-------------------------"
        cd tucker
        python -W ignore testTuckerDecomposition_withGreedy.py $2 $5
    fi
else
    help
    exit
fi
