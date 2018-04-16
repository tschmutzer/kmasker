# Installation
This instructions will help you to install kmasker including all prerequisites.
## Prerequisites
  To install and use kmasker you will need to following tools and libraries on your system:
- git
- blast
- jellyfish >= 2.2.6
- libjellyfish-devel
- libpthread (most likely already installed)
- [fastq-stats from ea-utils](https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/Compiling.md) 
- libgsl(-devel)
- [gffread](https://github.com/gpertea/gffread) 
- GCLIB as part of gffreads install instructions
- perl5
- scons
- tools for compiling c++ programs 
(we will try to offer binaries for differet OS to you)
## Installation Instructions
### Distribution specific parts
#### Ubuntu/Debian
 Installation of necessary tools and libs with apt
```bash
apt-get install jellyfish libjellyfish-2.0-dev gsl libgsl-dev scons blast2 build-essential
```

#### Redhat/CentOS
Installation of necessary tools and libs with yum
To be done. 

###General instructions (executed after the distrubution dependencies)
Create a directory for the source code of the different tools
`mkdir ~/src`
Check out the tools
```bash
cd ~/src
git clone https://github.com/ExpressionAnalysis/ea-utils
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread
git clone https://github.com/tschmutzer/kmasker
cd ea-utils/clipper
make fastq-stats
#Copy the binary file to a location which is in your path
#or to KMASKERs bin directory
cd ../..
cd gffread
make
cd ..
cd kmasker/src
scons
cp cmakser ../bin
cd ..
cp ../ea-utils/clipper/fastq-stats bin/
cp ../gffread/gffread bin/
#if you want you can add KMASKERs bin directory
#to your path
export PATH=$PATH:/$HOME/src/kmasker/bin
```
### MacOS
On MacOS the instructions are slightly different.
To be done.