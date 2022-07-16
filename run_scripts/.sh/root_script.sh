#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

# the directory of this run
run_dir=$lgame_home_dir/output/$run_id

if [ -d $run_dir ]
then
  rm -r $run_dir
fi
mkdir $run_dir

mv namelist.nml $run_dir

cd $run_dir

if [ ! -f ../../build/lgame ]
then
  echo "Executable lgame missing. Compile first. Aborting run."
  cd - > /dev/null
  exit 1
fi

if [ -f lgame ]
then
  rm lgame
fi

cp ../../build/lgame .

./lgame

cd - > /dev/null
