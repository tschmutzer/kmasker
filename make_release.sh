#!/bin/bash
git submodule init
git submodule update
./make_dockerimage.sh
./make_binaries.sh
./remove_dockerimage.sh
./make_structure.sh
