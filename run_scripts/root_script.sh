#!/bin/bash

# the directory of this run
run_dir=$rfpet_home_dir/output/$run_id

# creating the run directory if it does not exist
if [ ! -d $run_dir ]
then
mkdir $run_id
fi

cd $run_dir

cp ../bin/rfpet .
./rfpet


