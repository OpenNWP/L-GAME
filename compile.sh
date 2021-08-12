#!/bin/bash

# This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/RFPET

if [ ! -d build ]
then
mkdir build
fi

cd build

cmake -DCMAKE_INSTALL_PREFIX=../output ..
make

cd ..
