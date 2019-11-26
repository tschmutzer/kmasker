# Installation

These instructions will help you to install kmasker including all prerequisites.
You can install kmasker in the following ways:

- with bioconda in Linux or macOS (the easiest way)
- with precompiled binaries in Linux (easy)
- with the compilation of the source code in
    - Linux
    - macOS with the help of brew

# Installation with bioconda

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/kmasker/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/kmasker/badges/version.svg)](https://anaconda.org/bioconda/kmasker)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/kmasker/badges/latest_release_date.svg)](https://anaconda.org/bioconda/kmasker)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/kmasker/badges/platforms.svg)](https://anaconda.org/bioconda/kmasker)


*You should have installed bioconda before: [https://bioconda.github.io/](https://bioconda.github.io/)*
We recommend to create a new conda environment to avoid dependency issues. 
```bash
conda activate #if not active
conda create -n Kmasker Kmasker #create a environment called Kmasker
conda activate Kmasker #to use the Kmasker environment
#Check Kmasker installation
Kmasker --check_config --verbose
```

You can start kmasker with `Kmasker --help`
# Installation with precompiled binaries (Linux only)
We are working on an all-in-one Linux package for you,  with respect to the licences of the different tools. 
At the moment almost all tools are covered, except R. 
The binaries are located in the bin folder, their sources are located in the "external_src" folder.
You can find the releases at: 
[https://github.com/tschmutzer/kmasker/releases
](https://github.com/tschmutzer/kmasker/releases)
## 1. Distribution-depended instructions
The version with binaries has just a few requirements. Please consider the "Installation without prebuild binaries / newest version from git" if something is not working for you.

- Linux distribution not older than approx. 10 years 
- perl5
- R (tested with 3.5.0), but older/newer versions will work most likely
- Python (3.5) incl. pandas, numpy (if you want to use KRISPR module)

### Ubuntu/Debian
Installation of necessary tools and libs with apt 

```bash
sudo apt-get install r-base r-recommended 
```
If you want to use the KRISPR module you need the following python dependencies:

```bash
apt-get install python3.5 python3-pandas python3-numpy python3-pip
pip install kPAL
```
### Redhat/CentOS
Installation of necessary tools and libs with yum 

```bash
sudo yum install epel-release
sudo yum install R
```
If you want to use the KRISPR module please follow the RedHat developer instructions:
[https://developers.redhat.com/blog/2018/08/13/install-python3-rhel/
](https://developers.redhat.com/blog/2018/08/13/install-python3-rhel/)
and afterwards, run

```bash
pip install kPAL
```
## 2. General instructions (executed after the distribution dependencies)
Unzip the kmasker release.

```bash
unzip kmasker*.zip
cd kmasker
cd setup
./install_packages.R
#Check that all R packages were installed correctly
cd ..
#Optional: If you want to install kmasker, e.g. in /usr/local run
#scons PREFIX=/usr/local --release install
#(of course you need to install scons before)
```

You are done! Kmasker can be used now.

# Installation from source / newest version from git
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
- perl5 (with dependency data::Uniqid)
- scons
- tools for compiling c++ programs 
- R (tested with 3.5.0), but older/newer versions will work most likely
- Python (3.5) incl. pandas, numpy, kPAL (if you want to use KRISPR module)

## Installation Instructions
### 1. Distribution specific parts
#### Ubuntu/Debian
 Installation of necessary tools and libs with apt

```bash
apt-get install jellyfish libjellyfish-2.0-dev gsl libgsl-dev scons blast2 build-essential r-base r-recommended 
```
If you want to use the KRISPR module you need the following python dependencies:

```bash
apt-get install python3.5 python3-pandas python3-numpy python3-pip
pip install kPAL
```

#### Redhat/CentOS
##### Installation of necessary tools and libs with yum
Install development tools:

```bash
sudo yum group install "Development Tools"
```

Afterwards please compile Jellyfish from source: [https://github.com/gmarcais/Jellyfish]()
or download a binary release: [https://github.com/gmarcais/Jellyfish/releases]()

##### Install runtime dependencies:

```bash
sudo yum install epel-release
sudo yum install R
```
If you want to use the KRISPR module please follow the RedHat developer instructions:
[https://developers.redhat.com/blog/2018/08/13/install-python3-rhel/
]()
and afterwards, run

```bash
pip install kPAL
```

### 2. General instructions (executed after the distribution dependencies)
Create a directory for the source code of the different tools. 
For example in your home directory. 
`mkdir ~/src`
#### Check out the tools
```bash
cd ~/src
git clone https://github.com/ExpressionAnalysis/ea-utils
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread
git clone https://github.com/tschmutzer/kmasker #or unzip downloaded source code release!
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
#Download the submodule krispr
cd kmasker
git submodule init 
git submodule update
#Make cmasker
```
You can just compile cmasker and copy it to bin OR you can install Kmasker to a selected prefix.
If you have your jellyfish installation in a non-standard path, you can set JLIB and JINCLUDE.
For example `scons JLIB=/opt/jellyfish/lib/ JINCLUDE=/opt/jellyfish/includes/jellyfish-2.1.0/ -Q build`.
scons will also search in the selected prefix for jellyfish.
##### Just compile
```bash
scons -Q build
cp src/cmasker bin/cmasker
#if you want you, can add KMASKERs bin directory
#to your path
export PATH=$PATH:/$HOME/src/kmasker/bin
```
##### OR Install with a prefix
We can also compile and install Kmasker for you in a selected prefix, e.g. /usr/local.

```bash
scons PREFIX=/usr/local -Q install
#of course, you can also use JLIB and JINCLUDE here
#make sure that the PREFIX is in your PATH variable
```
#### Install R packages
```bash
cd setup
./install_packages.R
```


# macOS

**Information: You should be able to use the terminal in OSX. Please note that whether kmasker nor its dependencies has a graphical user interface.** 

Install brew from [https://brew.sh/](https://brew.sh/) (you can also use macports or fink, but this unsupported) and its dependencies (like OSX command line tools). Brew will guide you through its installation. We are working on a ready-to-use formula for brew at the moment. You can also use bioconda (instructions at the end of this document).
### Install Homebrew
Add the bioinformatics repository to brew.

```bash
brew tap brewsci/bio
brew install jellyfish blast scons
```
Download and install R from [https://cloud.r-project.org/bin/macosx/](https://cloud.r-project.org/bin/macosx/) or use brew cask:

```bash
brew install r-app
```

### Command line instructions
Create a directory for the source code of the different tools. 
For example in your home directory. 
`mkdir ~/src`
#### Check out the tools
```bash
cd ~/src
git clone https://github.com/ExpressionAnalysis/ea-utils
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread
git clone https://github.com/tschmutzer/kmasker #or unzip downloaded source code release!
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
#Download the submodule krispr
cd kmasker
git submodule init 
git submodule update
#Make cmasker
#You can get the current version with brew info jellyfish
scons JINCLUDE=/usr/local/include/jellyfish-2.2.10 JLIB=/usr/local/lib -Q build
cp src/cmasker bin/cmasker
#or install with 
#scons PREFIX=/your/prefix JINCLUDE=/usr/local/include/jellyfish-2.2.10 JLIB=/usr/local/lib  -Q install
```
#### Install R packages
```bash
cd setup
./install_packages.R
```

#### Install python3 dependencies (KRISPR)
```bash
brew install python3
pip3 install numpy
pip3 install pandas
pip3 install kPAL

```

You are done. You can use Kmasker with ~/src/kmasker/bin/Kmasker or in your prefix. 

