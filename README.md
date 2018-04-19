# Kmasker

Usage of program Kmasker:
 (version:  0.0.24 rc180417)

## Description:

Kmasker is a tool for the automatic detection of repetitive sequence regions.

### Modules:

 --build                 construction of new index (requires --indexfiles)
 
 --run                   run k-mer repeat detection and masking (requires --fasta)
 
 --explore               perform downstream analysis with constructed index and detected repeats
 

### General options:

 --show_repository       shows complete list of global and private k-mer indices
 
 --show_details          shows details for a requested kindex
 
 --remove_kindex         remove kindex from repository
 
 --expert_setting        submit individual parameter to Kmasker (e.g. on memory usage for index construction)


### Installation & Requirements:

See our Wiki for details on how to install Kmasker. There, you find our list of requirements of external tools. 
Please make sure, that these are in your PATH environemnt. If not please specify them in the 'kmasker.config' file.

Kmasker uses in interal repository for reuse of kindex structures. The data (calculated kindex) will be stored either in lokal ('private') or in the global directory. The path to the global directory has to be set in the kmasker.config file after installation.


## Commands:

### Quick command overview:
Kmasker --help

Kmasker --build --seq sequence.fastq --gs 135 --in At1 --cn arabidopsis

Kmasker --run --fasta query.fasta --kindex At1

Kmasker --explore --hexplot --multi_kindex At1 Hv1

Kmasker --show_repository

Kmasker --show_details At1


### [BUILD]:

The build module is used to construct a k-mer index structure. It has its own help section (type '--help'). 
One either can provide parameters using the command line or use the option '--config' to provide a config file with detailed meta data. If parameters are defined twice, the config file will overwrite parameters given at the command line.

### [RUN]:

The run module starts the core process of Kmasker. There are two general options. 1.) Analyse input with SINGLE k-mer index structures and 2) perform a comparative analysis using MULTIPLE (2) k-mer index structures.

#### SINGLE

Kmasker --run --fasta query.fasta --kindex At1

#### MULTIPLE

Kmasker --run --fasta query.fasta --multi_kindex At1 Hv1

### [EXPLORE]:

The explore module provides additional functionality for downstream analysisi e.g. vizualisations or annotation. 

#### EXPLORE - ANNOTATION

Kmasker --explore --annotate --fasta query.fasta --gff kmasker_result.gff --dbfasta mipsREdat_9.3p_ALL.fasta

For repeat studies in plant species we recommend using curated libraries:

PGSB-REdat
Link: http://pgsb.helmholtz-muenchen.de/plant/recat/

TREP:
http://botserv2.uzh.ch/kelldata/trep-db/


#### EXPLORE - VISUALISATION


HEXPLOT

Kmasker --explore --hexplot --multi_kindex At1 Hv1

This visualisation is can be used for comparitve inspection of 2(!) k-mer index structures (e.g. to study two species by comparing the correpsonding WGS data).


TRIPLOT

Kmasker --explore --triplot --multi_kindex At1 Hv1 Sc1

This visualisation is comparing 3(!) k-mer index structures in a triangular plot using the R package ('ggtern' [LINK]). It provides a broader perspective which k-mers are shared with equal (centered) or different (outer regions) abundances throuout the investigated datasets (e.g. three different species).


