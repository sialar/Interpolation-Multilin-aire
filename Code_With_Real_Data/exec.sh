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
#"$PROJECT_PATH"/bin/TestAdaptativeInterpolation $*

core=$2

if [ $# -le 2 ]
then
    exit
fi

shift
shift

cd AI/python

if [ $1 = "ALL" ]
then
    args="macro_totale0 macro_totale1 macro_fission0 macro_fission1 macro_nu_fission0 "
    args="$args macro_nu_fission1 macro_scattering000 macro_scattering001 macro_scattering010 "
    args="$args macro_scattering011 macro_absorption0 macro_absorption1"
else
    args=$*
fi

for arg in $args
do
    python plot_results.py $core $arg
done

if [ $1 = "ALL" ]
then
  python plot_reactivity.py $core
fi

exit
