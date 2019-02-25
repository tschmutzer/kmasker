# Kmasker

Usage of program Kmasker:
 (version:  0.0.35 rc190108)

## Description:

Kmasker is a tool for the automatic detection of repetitive sequence regions.

### Modules:

 --build                 construction of new index
 
 --run                   perform k-mer repeat detection, masking, abundance and comparative analysis
 
 --explore               perform downstream analysis with constructed index and detected repeats
 

### General options:

 --show_repository       shows complete list of global and private k-mer indices
 
 --show_details          shows details for a requested kindex
 
 --show_path             shows path like private and external path which are checked for pre-computed kindex structure
 
 --remove_kindex         remove kindex from repository
 
 --expert_setting_*      submit individual parameters to subprocesses (* can be blast, jellyfish or kmasker).
 
 --config_*              submit individual parameters in a configuration file (* can be blast, jellyfish or kmasker).


### Installation & Requirements:

See our [Wiki](https://github.com/tschmutzer/kmasker/wiki/01-Installation) for details on how to install Kmasker. There, you find our list of requirements of external tools. 
Please make sure, that these are in your PATH environemnt. If not please specify them in the 'kmasker.config' file.

Kmasker uses an interal repository for reuse of kindex structures. The data (calculated kindex) will be stored either in local ('private') or in the external directory. The path to the external directory has to be set in the kmasker.config file after installation.


## Commands:

### Quick command overview:

Kmasker --help

Kmasker --show_repository

Kmasker --show_details At1

Kmasker --build --seq sequence.fastq --gs 135 --in At1

Kmasker --run --fasta query.fasta --kindex At1

Kmasker --explore --hexplot --kindex At1 Hv1



## [BUILD]:

The build module is used to construct a k-mer index structure. It has its own help section (type '--help'). 
One either can provide parameters using the command line or use the option '--config' to provide a config file with detailed meta data. If parameters are defined twice, the config file will overwrite parameters given at the command line.

Kmasker --build --seq sequence.fastq --gs 135 --in At1 --cn arabidopsis



## [RUN]:

The run module starts the core process of Kmasker. There are four general options: 1) basic k-mer analysis with SINGLE or MULTIPLE index structures, 2) screen a sequence for candidates applicable in 'fluorescent in situ hybridization' (FISH) 3) perform a comparative analysis searching for differences in applied k-mer index structures and 4) analysis of short sequence probes for its genome-wide specificity.

#### SINGLE or MULTIPLE

Kmasker --run --fasta query.fasta --kindex At1

Kmasker --run --fasta query.fasta --kindex At1 Hv1

#### FISH

Kmasker --run --fish --fasta query.fasta --kindex At1

Check input sequence for long sequence stretches with low repetitiveness. As a results candidate sequences with good target specificity and functionality are selected.

#### KRISPR

Kmasker --run --krispr --fasta candidate_grna.fasta --kindex Hv1

Check a set of gRNAs for their specificity in CRISPR/cas application. As a result for each given gRNA a score is calculated reflecting its target specififity. 

#### COMPARATIVE

Kmasker --run --compare --fasta query.fasta --kindex At1 Hv1

Perform a comparative study using multiple kindex. Kmasker will search for differences between both applied kindex.


## [EXPLORE]:

The explore module provides additional functionality for downstream analysis e.g. vizualisations or annotation. 

#### EXPLORE - VISUALISATION


HISTPLOT

Two kinds of histograms can be constructed ('hist' and 'histm'). 

Kmasker --explore --hist --occ species_A.occ --list selection.ids

Kmasker --explore --histm --occ species_A.occ --list selection.ids

The first histogram type ('hist') will construct a graphic that shows raw k-mer frequencies per sequence. The second histogram type ('histm') will construct a graphic with mean values calculated in a sliding window approach. Furthermore, the 'histm' plots use log scales for improved visability. If large dataset with millions of contigs are analyzed its recommended to use both methods with the '--list' parameter that is providing as subset of selected sequence identifiers. This will avoid long computing times. 


HEXPLOT

Kmasker --explore --hexplot --occ species_A.occ species_B.occ

Kmasker --explore --hexplot --file comparative_analysis.stats

This visualisation can be used for comparative inspection of two constructed Kmasker output files (occ). Keep in mind that only sequences are compared that are present in both files. The visualisation illustrates which sequences are the most different based on k-mer frequencies (e.g. to study two species by comparing the correpsonding WGS data).



#### EXPLORE - ANNOTATION

Kmasker --explore --annotate --fasta query.fasta --gff kmasker_result.gff --feature KRC --dbfasta mipsREdat_9.3p_ALL.fasta

For repeat studies in plant species we recommend using curated libraries:

PGSB-REdat
Link: http://pgsb.helmholtz-muenchen.de/plant/recat/

TREP:
http://botserv2.uzh.ch/kelldata/trep-db/


