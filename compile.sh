#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

if [ ! -d build ]
then
  mkdir build
fi

cd build

d_value=False
f_value=False
while getopts "df" opt; do
  case $opt in
    d)
      d_value=True
      ;;
    f)
      f_value=True
      ;;
    \?)
      echo "Invalid option: -$OPTARG. Compiling anyway."
      ;;
  esac
done

cmake -DBOUNDS_CHECKS=$d_value -DFAST=$f_value ..
make

cd ..
