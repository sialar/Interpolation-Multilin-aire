#!/bin/sh


PROJECT_PATH="/home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data/AI"
make -j8

help()
{
  echo ""
  echo " Arguments required:"
  echo "   - arg 1  : Number of iteration in AI algorithm [$1]"
  echo "   - arg 2  : Core type [$2]"
  echo "   - arg >2 : Reaction types (between 1 and 12) [$2]"
  echo ""
}

help

if [ "$1" = "AI" ]
then
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

elif [ "$1" = "Convergence" ]
then
    shift
    "$PROJECT_PATH"/bin/TestConvergence $*

    cd AI/python
    python plot_errors.py $*
fi

exit
