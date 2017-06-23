#!/bin/sh


PROJECT_PATH="/home/sialar/Stage/LaboJ_LLions/Code/Code_EDF_Data/AI"
make -j8

help()
{
  echo ""
  echo " To launch the script, you need to choose one argument from the list below:"
}

help

"$PROJECT_PATH"/bin/TestAdaptativeInterpolation $1

exit
