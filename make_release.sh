#!/bin/bash
./make_dockerimage.sh
./make_binaries.sh
./remove_dockerimage.sh
./make_structure.sh
