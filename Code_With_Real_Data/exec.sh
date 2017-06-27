#!/bin/sh


PROJECT_PATH="/home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data/AI"
make -j8

help()
{
  echo ""
  echo " Arguments required:"
  echo "   - arg 1  : Number of iteration in AI algorithm [$1]"
  echo "   - arg 2  : Core type [$2]"
  echo "   - arg >2 : reaction types (between 1 and 12) [$2]"
  echo ""
}

help

"$PROJECT_PATH"/bin/TestAdaptativeInterpolation $*

exit
