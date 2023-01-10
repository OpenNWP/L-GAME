#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

# downloading GLDB data

wget "http://www.flake.igb-berlin.de/data/gldbv2.tar.gz"
tar -xzf gldbv2.tar.gz GlobalLakeDepth.dat
rm gldbv2.tar.gz
