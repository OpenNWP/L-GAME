#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

if [ ! -d output ]
then
mkdir output
fi

if [ ! -d real_weather ]
then
mkdir real_weather
fi

if [ ! -d figs ]
then
mkdir figs
fi

if [ ! -d grids ]
then
mkdir grids
fi

if [ ! -d grids/phys_sfc_properties ]
then
mkdir grids/phys_sfc_properties
fi
