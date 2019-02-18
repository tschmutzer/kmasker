#!/bin/bash

set -e

source /hbb_exe/activate
cd /io/external/jellyfish/
env CFLAGS="$STATICLIB_CFLAGS" CXXFLAGS="$STATICLIB_CXXFLAGS"   ./configure --prefix=/hbb_exe/ --disable-shared --enable-static
make
make install
#ln -s /hbb_exe/include/jellyfish-2.2.10/jellyfish /hbb_exe/include/jellyfish

cd /io/src
scons

exit
