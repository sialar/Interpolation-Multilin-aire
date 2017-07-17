#!/bin/sh


PROJECT_PATH="/home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data_N1/AI"
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

NITER=$1
CORE=$2

if [ $# -le 2 ]
then
    exit
fi

shift
shift

if [ $1 = "ALL" ]
then
    ARGS="macro_totale0 macro_totale1 macro_fission0 macro_fission1 macro_nu*fission0 "
    ARGS="$ARGS macro_nu*fission1 macro_scattering000 macro_scattering001 macro_scattering010 "
    ARGS="$ARGS macro_scattering011 macro_absorption0 macro_absorption1"
else
    ARGS=$*
fi

for arg in $ARGS
do
    "$PROJECT_PATH"/bin/TestAdaptativeInterpolation $NITER $CORE $arg
done

cd AI/python

for arg in $ARGS
do
    python plot_results.py $CORE $arg
done

REACTIVITY_CS="macro_totale0 macro_totale1 macro_nu_fission0 macro_nu_fission1 "
REACTIVITY_CS="$REACTIVITY_CS macro_scattering000 macro_scattering001 "
REACTIVITY_CS="$REACTIVITY_CS macro_scattering010 macro_scattering011"

for cs in $REACTIVITY_CS
do
    file="$PROJECT_PATH/data/$CORE/$cs"
    if [ ! -f "$file" ]
    then
        exit
    fi
done
#python plot_reactivity.py $CORE

exit
