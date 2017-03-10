# Kmasker

Usage of program Kmasker:
 (version:  0.0.21 rc170114)

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

