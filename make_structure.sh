#!/bin/bash
mv src/cmasker bin/cmasker
mv external_src/fastq-stats bin/
mv external_src/jellyfish/bin/jellyfish bin/
mv external_src/gffread_bin bin/gffread
curl "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$(curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ --list-only | grep ncbi-blast-.\*-linux\\.tar\\.gz\$)" -o blast.tar.gz
tar -xzvf blast.tar.gz ncbi-blast\*/bin/makeblastdb ncbi-blast\*/bin/blastn
cp ncbi-blast*/bin/makeblastdb bin/
cp ncbi-blast*/bin/blastn bin/
rm -rf ncbi-blast*
rm -rf blast.tar.gz
cd ..
zip -9 -r --exclude=\*.git\* --exclude=\*.DS_Store\*  --exclude=\*make_\*.sh --exclude=\*remove_dockerimage.sh --exclude=\*/builder/\* --exclude=\*/example/\*  kmasker_release.zip kmasker
