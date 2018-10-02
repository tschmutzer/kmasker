#!/bin/bash
cp -R external_src/jellyfish-2.2.10 external_src/jellyfish
docker run -t -i --rm -v `pwd`/src:/io/src -v `pwd`/external_src:/io/external kmasker/builder bash /io/external/make_jellyfish_cmakser.sh
docker run -t -i --rm -v `pwd`/src:/io/src -v `pwd`/external_src:/io/external kmasker/builder bash /io/external/make_fastqstats.sh
docker run -t -i --rm -v `pwd`/src:/io/src -v `pwd`/external_src:/io/external kmasker/builder bash /io/external/make_gffread.sh
