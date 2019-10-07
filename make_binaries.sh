#!/bin/bash
wget -c https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
tar -xvf jellyfish-2.2.10.tar.gz
mv jellyfish-2.2.10 external_src/jellyfish
rm -rf jellyfish-2.2.10.tar.gz
docker run -t -i --rm -v `pwd`/:/io/src -v `pwd`/external_src:/io/external kmasker/builder bash /io/external/make_jellyfish_cmakser.sh
docker run -t -i --rm -v `pwd`/external_src:/io/external kmasker/builder bash /io/external/make_fastqstats.sh
docker run -t -i --rm -v `pwd`/external_src:/io/external kmasker/builder bash /io/external/make_gffread.sh
