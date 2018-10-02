#!/bin/bash
cd /io/external/gffread
source /hbb_exe/activate
make
cd /io/external
cp ./gffread/gffread gffread_bin
cd gffread
make clean
exit
