# Installation
This instructions will help you to install kmasker including all prerequisites.
# Installtion from release with binaries (linux only)
We are working on an all-in-one linux package for you taking account the licences of the different tools. 
At the moment almost all tools are covered, except R (we are working on that in the moment). 
The binaries are located in the bin folder, their sources are located in the "external_src" folder.
## Prerequisites
The version with binaries has just a few requirements. Please consinder the "Installation without prebuild binaries / newest version from git" if something is not working for you.

- Linux distribution not older than approx. 10 years 
- perl5
- R (tested with 3.5.0), but older/newer versions will work most likely

#### Ubuntu/Debian
 Installation of necessary tools and libs with apt
```bash
sudo apt-get install r-base r-recommended 
```
#### Redhat/CentOS
```bash
sudo yum install epel-release
sudo yum install R
```

### General instructions (executed after the distrubution dependencies)
Unzip the kmasker release.
```bash
unzip kmasker*.zip
cd kmasker
cd setup
./install_packages.R
```

Your are done! Kmasker can be used now.

# Installation without prebuild binaries / newest version from git / Macintosh
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
- R (tested with 3.5.0), but older/newer versions will work most likely
(we will try to offer binaries for different OS to you)
## Installation Instructions
### Distribution specific parts
#### Ubuntu/Debian
 Installation of necessary tools and libs with apt
```bash
apt-get install jellyfish libjellyfish-2.0-dev gsl libgsl-dev scons blast2 build-essential r-base r-recommended 
```

#### Redhat/CentOS
Installation of necessary tools and libs with yum
To be done. 

### General instructions (executed after the distrubution dependencies)
Create a directory for the source code of the different tools. 
For example in our home directory. 
`mkdir ~/src`
#### Check out the tools
```bash
cd ~/src
git clone https://github.com/ExpressionAnalysis/ea-utils
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread
git clone https://github.com/tschmutzer/kmasker
```
#### Compile the tools
```bash
cd ea-utils/clipper
make fastq-stats
#Copy the binary file to a location which is in your path
#or to KMASKERs bin directory
cp fastq-stats ../../kmasker/bin/
cd ../..
cd gffread
make
cp gffread ../kmasker/bin/
cd ..
cd kmasker/src
scons
cp cmasker ../bin
cd ..
#if you want you can add KMASKERs bin directory
#to your path
export PATH=$PATH:/$HOME/src/kmasker/bin
```
#### Install R packages
```bash
cd setup
./install_packages.R
```

### MacOS

**Information: You should be able to use the terminal in OSX. Please note that whether kmasker nor its depedencies has a graphical user interface.** 

Install brew from https://brew.sh/ (you can also use macports or fink, but this unsupported) and its depenedencies (like OSX command line tools). Brew will guide you through its installation. 
#### Install prerequisites
Add the bioinformatics repository to brew.
```bash
brew tap brewsci/bio
brew install jellyfish blast scons
```
Download and install R from https://cloud.r-project.org/bin/macosx/

### General instructions (executed after the distrubution dependencies)
Create a directory for the source code of the different tools. 
For example in our home directory. 
`mkdir ~/src`
#### Check out the tools
```bash
cd ~/src
git clone https://github.com/ExpressionAnalysis/ea-utils
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread
git clone https://github.com/tschmutzer/kmasker
```
#### Compile the tools
```bash
cd ea-utils/clipper
make fastq-stats
#Copy the binary file to a location which is in your path
#or to KMASKERs bin directory
cp fastq-stats ../../kmasker/bin/
cd ../..
cd gffread
make
cp gffread ../kmasker/bin/
cd ..
cd kmasker/src
export CXXFLAGS="-I/usr/local/include/jellyfish-2.2.10 -L/usr/local/lib $CXXFLAGS"
#Modify the above command to fit your jellyfish version
#You can get the current version with brew info jellyfish
scons
cp cmasker ../bin
cd..
```
#### Install R packages
```bash
cd setup
./install_packages.R
```
You are done. You can use Kmasker with ~/src/kmasker/bin/Kmasker
