# Kmasker

Usage of program Kmasker:
 (version:  0.0.25 rc180723)

## Description:

Kmasker is a tool for the automatic detection of repetitive sequence regions.

### Modules:

 --build                 construction of new index
 
 --run                   run k-mer repeat detection, masking, abundance and comparative analysis
 
 --explore               perform downstream analysis with constructed index and detected repeats
 

### General options:

 --show_repository       shows complete list of global and private k-mer indices
 
 --show_details          shows details for a requested kindex
 
 --show_path             shows path like private and external path which are checked for pre-computed kindex structure
 
 --remove_kindex         remove kindex from repository
 
 --expert_setting_*      submit individual parameters to subprocesses (* can be blast, jellyfish or kmasker).
 
 --config_*              submit individual parameters in a configuration file (* can be blast, jellyfish or kmasker).


### Installation & Requirements:

See our Wiki for details on how to install Kmasker. There, you find our list of requirements of external tools. 
Please make sure, that these are in your PATH environemnt. If not please specify them in the 'kmasker.config' file.

Kmasker uses an interal repository for reuse of kindex structures. The data (calculated kindex) will be stored either in local ('private') or in the external directory. The path to the external directory has to be set in the kmasker.config file after installation.


## Commands:

### Quick command overview:

Kmasker --help

Kmasker --show_repository

Kmasker --show_details At1

Kmasker --build --seq sequence.fastq --gs 135 --in At1

Kmasker --run --fasta query.fasta --kindex At1

Kmasker --explore --hexplot --multi_kindex At1 Hv1



## [BUILD]:

The build module is used to construct a k-mer index structure. It has its own help section (type '--help'). 
One either can provide parameters using the command line or use the option '--config' to provide a config file with detailed meta data. If parameters are defined twice, the config file will overwrite parameters given at the command line.

Kmasker --build --seq sequence.fastq --gs 135 --in At1 --cn arabidopsis



## [RUN]:

The run module starts the core process of Kmasker. There are three general options. 1.) Analyse input with SINGLE k-mer index structures, 2) perform a comparative analysis using MULTIPLE (2) k-mer index structures and 3) analysis of short sequence probes for its genome-wide specificity.

#### SINGLE

Kmasker --run --fasta query.fasta --kindex At1

#### MULTIPLE

Kmasker --run --fasta query.fasta --multi_kindex At1 Hv1

#### gRNA

Kmasker --run --grna candidate_gRNA.fasta --kindex Hv1

Check a set of gRNAs for their specificity in CRISPR/cas application. As a result for each given gRNA a score is calculated reflecting its target specififity. 


## [EXPLORE]:

The explore module provides additional functionality for downstream analysisi e.g. vizualisations or annotation. 

#### EXPLORE - VISUALISATION


HISTPLOT

Two kinds of histograms can be constructed ('hist' and 'histm'). 

Kmasker --explore --hist --occ species_A.occ --list selection.ids

Kmasker --explore --histm --occ species_A.occ --list selection.ids

The first type ('hist') will construct histograms showing raw k-mer frequencies and calculated means etc. per sequence. The second type ('histm') will construct histograms using mean values calculated in a sliding window approach. Furthermore, these plots use log scales for improved visability. If large dataset with millions of contigs are analyzed its recommended to use both methods with the '--list' parameter that is providing as subset of selected sequence identifiers. This will avoid long computing times. 


HEXPLOT

Kmasker --explore --hexplot --occ species_A.occ species_B.occ

This visualisation can be used for comparative inspection of two constructed Kmasker output files (occ). Keep in mind that only sequences are compared that are present in both files. The visualisation illustrates which sequences are the most different based on k-mer frequencies (e.g. to study two species by comparing the correpsonding WGS data).



#### EXPLORE - ANNOTATION

Kmasker --explore --annotate --fasta query.fasta --gff kmasker_result.gff --feature KRC --dbfasta mipsREdat_9.3p_ALL.fasta

For repeat studies in plant species we recommend using curated libraries:

PGSB-REdat
Link: http://pgsb.helmholtz-muenchen.de/plant/recat/

TREP:
http://botserv2.uzh.ch/kelldata/trep-db/


