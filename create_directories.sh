#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

if [ ! -d real_weather ]
then
  mkdir real_weather
fi

if [ ! -d nwp_init ]
then
  mkdir nwp_init
fi
