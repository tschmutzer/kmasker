# Kmasker

Usage of program Kmasker:
 (version:  0.0.23 rc170504)

## Description:

Kmasker is a tool for the automatic detection of repetitive sequence regions.

### Modules:

 --build                 construction of new index (requires --indexfiles)
 --run                   run k-mer repeat detection and masking (requires --fasta)
 --postprocessing        perform downstream analysis with constructed index and detected repeats

### General options:

 --show_repository       shows complete list of global and private k-mer indices
 --show_details          shows details for a requested kindex
 --remove_kindex         remove kindex from repository
 --expert_setting        submit individual parameter to Kmasker (e.g. on memory usage for index construction)


### Requirements:
Here, we provide a list of external tools that are used within Kmasker. Please make sure, that these are in your PATH environemnt. If not please specify them in the 'config.kmasker' file.

## Commands:

### Quick command overview:
Kmasker --help
Kmasker --build --seq sequence.fastq --gs 135 --name At1
Kmasker --run --fasta query.fasta --kindex At1
Kmasker --show_repository
Kmasker --show_details At1


### --build:

The build module is used to construct a k-mer index structure. It has its own help section (type '--help'). 
One either can provide parameters using the command line or use the option '--config' to provide a config file with detailed meta data.

### --run:

The run module starts the core process of Kmasker. Your input sequence is provided 

There are two general options. 1.) Analyze input with SINGLE k-mer index structures and 2) perform a comparative analysis using MULTIPLE (2) k-mer index structures.

### --postprocessing:

The postprocessing module provides additional analysis that can be performed using the constructed KINDEX structures e.g. vizualisations.
