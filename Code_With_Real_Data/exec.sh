#!/bin/sh


PROJECT_PATH="/home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data/AI"
make -j8

helpAI()
{
  echo ""
  echo " Arguments required:"
  echo "   - arg 1  : Number of iteration in AI algorithm [$1]"
  echo "   - arg 2  : Core type [$2]"
  echo "   - arg >2 : Reaction types (between 1 and 12) [$2]"
  echo ""
}

helpConvergence()
{
  echo ""
  echo " Arguments required:"
  echo "   - arg 1  : Core type [$2]"
  echo "   - arg >1 : Reaction types (between 1 and 12) [$2]"
  echo ""
}

help()
{
  echo ""
  echo " Arguments required:"
  echo "   - arg 1 : --> AI to launch Adaptative Interpolation algorithm"
  echo "             --> Convergence to get the convergence speed graph"
  echo ""
}

if [ "$1" = "AI" ]
then
    helpAI
    shift
    "$PROJECT_PATH"/bin/TestAdaptativeInterpolation $*
    core=$2

    if [ $# -le 2 ]
    then
        exit
    fi

    shift
    shift

    cd AI/python

    if [ "$1" = "ALL" ]
    then
        args="macro_totale0 macro_totale1 macro_absorption0 macro_absorption1 "
        args="$args macro_scattering000 macro_scattering001 macro_scattering010 macro_scattering011"
        args="$args macro_nu_fission0 macro_nu_fission1 macro_fission0 macro_fission1"
    else
        args=$*
    fi

    for arg in $args
    do
        python plot_results.py $core $arg
        python plot_results_with_cocagne.py $core $arg
    done

    if [ "$1" = "ALL" ]
    then
        python plot_reactivity.py $core
        python plot_reactivity_with_cocagne.py $core
    fi
    #for i in `seq 0 4`;
    #do
    #    python plot_interpolation_points_1D.py $core $i
    #done

    #for i in `seq 0 4`;
    #do
    #    for j in `seq 0 4`;
    #    do
    #        if [ $i -lt $j ]
    #        then
    #            python plot_interpolation_points.py $core $i $j
    #        fi
    #    done
    #done


elif [ "$1" = "Convergence" ]
then
    helpConvergence
    shift
    "$PROJECT_PATH"/bin/TestConvergence $*

    cd AI/python
    python plot_errors.py $*

else
    help
    exit
fi
