#!/bin/bash

# the directory of this run
run_dir=$rfpet_home_dir/output/$run_id

# creating the run directory if it does not exist
if [ ! -d $run_dir ]
then
mkdir $run_dir
fi

mv namelist.nml $run_dir

cd $run_dir

if [ ! -f ../../build/rfpet ]
then
echo "Executable rfpet missing. Compile first. Aborting run."
cd - > /dev/null
exit 1
fi

if [ -f rfpet ]
then
rm rfpet
fi

cp ../../build/rfpet .

mpirun -np $ncpus ./rfpet

cd - > /dev/null
