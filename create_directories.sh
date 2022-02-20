#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

if [ ! -d output ]
then
mkdir output
fi

if [ ! -d bc ]
then
mkdir bc
fi

if [ ! -d figs ]
then
mkdir figs
fi