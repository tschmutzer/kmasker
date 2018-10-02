#!/bin/bash
cd /io/external/ea-utils/clipper
source /hbb_exe/activate
export CPLUS_INCLUDE_PATH=/io/external/ea-utils/clipper:/hbb_exe/include
export C_INCLUDE_PATH=/io/external/ea-utils/clipper:/hbb_exe/include
make fastq-stats
cp fastq-stats ../../fastq-stats
make clean
exit
