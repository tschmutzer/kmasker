# Kmasker

## Installation
The installation instructions are available at: [INSTALL.md](https://github.com/tschmutzer/kmasker/blob/master/INSTALL.md). There, you find our list of requirements of external tools. Please make sure, that these are in your PATH environemnt if you installed from source. If not please specify them in the 'etc/kmasker.config' file.

Kmasker uses an interal repository for reuse of kindex structures. The data (calculated kindex) will be stored either in local ('private' - definied in ~/.kmasker_user.config) or in the external directory. The path to the external directory has to be set in the 'etc/kmasker.config' file after installation.

## Description:

Kmasker plants is a tool for the automatic detection of sequence regions with meaningful k-mer characteristics. This can be sequences with highly abundant k-mer patterns (repeats), regions with diverging k-mer patterns between two studied WGS samples or segments with high target specificity.

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


## Config files
Kmasker has two config files. A global one which is in prefix/etc/kmasker.conf and a local one which is in ~/.kmasker_user.conf
Both files define the same parameters, but the local parameters will overrule the global ones. 
The global conf file is inted for standard parameters on a multi-user system (like a server). But sometimes it is useful to define parameters for yourself in the local file.
### Global
Lets have a look on the global config: 

```
#Kmasker requirements
PATH_kindex_global=/hsm/hsm/schmutzr/KINDEX/

#external tools
jellyfish=
fastq-stats=
gffread= 
```
You can define the system-wide settings here. PATH_kindex_global will store the large kindex (jellyfish indicies) files, thus it should be on a volume with plenty of hard disk space. 
If the external tools are not in the default path (for every user), they can be definied here. A `kmasker --check_install` will check them and also rewrite them based on your path if the entries are empty. This command should be executed by the administrator. 
### User
The user config can be basically the same. You can even override the global kindex path (e.g. if different indicies for different projects are provided). The biggest difference is the PATH\_kindex\_private path. Here you can provide a personal repository (e.g. for testing a new organism, genotype etc.) for building new KINDEX indicies. This is espacially useful if you are not allowed to write in the global path. 
Also you can overwrite the tool paths, e.g. if you want to use a different R version (don't forget to reinstall the packages).

```
PATH_kindex_private=/mount/myvolume/KINDEX
#external tool requirements
gffread=/Users/unibenutzer/src/gffread/gffread
R=/usr/local/bin/R
makeblastdb=/usr/local/bin/makeblastdb
jellyfish=/usr/local/bin/jellyfish
blastn=/usr/local/bin/blastn
fastq-stats=/Users/chrisulpinnis/src/ea-utils/clipper/fastq-stats
```

### Kmasker commands with config files
`--check_install` will check your installation and create the global config file (or fill missing values in it).

From the tool's help:
`--user_conf			 set specific user configuration file [/Users/chrisulpinnis/.kmasker_user.config]`
`--global_conf			 set specific global configuration file [/Users/chrisulpinnis/miniconda3/share/kmasker/etc/kmasker.config]`
`--show_path			 show path Kmaskers looks for constructed kindex`

## Commands:

### Quick command overview:

Kmasker --help

Kmasker --show_repository

Kmasker --show_details At1

Kmasker --build --seq sequence.fastq --gs 135 --in At1

Kmasker --run --fasta query.fasta --kindex At1

Kmasker --explore --hexplot --kindex At1 Hv1

The build in help can be used in general:

```
Kmasker --help

Usage of program Kmasker:
 (version:  0.0.36 rc191017)				(session id: 8sAmvoVLq)

 Description:
	 Kmasker is a tool for the automatic detection of repetitive sequence regions.
	 There are three modules and you should select one for your analysis.

 Modules:
 --build		 construction of new index (requires --seq)
 --run			 perform analysis and masking (requires --fasta)
 --explore		 perform downstream analysis with constructed index and detected repeats

 General options:
 --show_repository		 show complete list of private and external k-mer indices
 --show_details			 show details for a requested kindex
 --show_path			 show path Kmaskers looks for constructed kindex
 --remove_kindex		 remove kindex from repository
 --set_private_path		 change path to private repository
 --set_external_path		 change path to external repository [readonly]
 --expert_setting_kmasker	 submit individual parameter to Kmasker eg. pctgap,
				 minseed, mingff (see documentation!)
 --expert_setting_jelly		 submit individual parameter to jellyfish (e.g. on memory usage
				 for index construction)
 --expert_setting_blast		 submit individual parameter to blast (e.g. '-evalue')
 --threads			 set number of threads [4]
 --bed				 force additional BED output [off]
 --user_conf			 set specific user configuration file [/Users/unibenutzer/.kmasker_user.config]
 --global_conf			 set specific global configuration file [/Users/unibenutzer/miniconda3/envs/kmasker/etc/kmasker.config]
 --check_install		 shows the detected/configured path for all used applications
 --setid			 set a user specified process id
 --long_id			 create a process id that is unique for this host (e.g. for use in cluster environments)
 --temp				 sets the location of temporary files [./temp/]
 --verbose			 enables verbose output and keeps log files
```
or for each module.



## [BUILD]:

The build module is used to construct a k-mer index structure. It has its own help section (type '--help'). 
One either can provide parameters using the command line or use the option '--config' to provide a config file with detailed meta data. If parameters are defined twice, the config file will overwrite parameters given at the command line.

Kmasker --build --seq sequence.fastq --gs 135 --in At1 --cn arabidopsis

### Expert settings
Sometimes a server or desktop computer has not enough computational power to calculate the indicies with jellyfish. 

E.g. for a desktop computer with 12 threads and 16GB memory:

`Kmasker --build --seq assembly3_WGSMorex_renamed_blastable_carma.fasta -gs 5 --in WGSMorex --cn barley --verbose --threads 12 --expert_setting_jelly size=3`
The `--expert_setting_jelly size=3` will add --size 3 to jellyfish. You can also specifiy the threads (regardless of Kmaskers --threads flag - in case you want to change the number of threads just for jellyfish).
### Module help
```
Kmasker --help --build

Usage of program Kmasker:
 (version:  0.0.36 rc191017)				(session id: 4rrLaaftil)

 Command:
	 Kmasker --build --seq mysequences.fasta

 Options:
 --seq		 fasta or fastq sequence(s) that are used to build the index
 --k		 k-mer size to build index [21]
 --gs		 genome size of species (in Mbp)
 --in 		 provide k-mer index name (e.g. HvMRX for hordeum vulgare cultivare morex) [date]
 --cn 		 provide common name of species (e.g. barley)
 --as 		 input is of sequence type 'assembly'. Set this option if your input is based on
      		 assembled sequences (e.g. genome reference). Default sequence type is 'reads'.
```



## [RUN]:

The run module starts the core process of Kmasker. There are four general options: 1) basic k-mer analysis with SINGLE or MULTIPLE index structures, 2) screen a sequence for candidates applicable in 'fluorescent in situ hybridization' (FISH) 3) perform a comparative analysis searching for differences in applied k-mer index structures and 4) analysis of short sequence probes for its genome-wide specificity.

#### SINGLE or MULTIPLE

`Kmasker --run --fasta query.fasta --kindex At1`

`Kmasker --run --fasta query.fasta --kindex At1 Hv1`

#### FISH

`Kmasker --run --fish --fasta query.fasta --kindex At1`

Check input sequence for long sequence stretches with low repetitiveness. As a results candidate sequences with good target specificity and functionality are selected.

#### KRISPR

##### Use krispr with default model

`Kmasker --run --krispr --fasta candidate_grna.fasta --kindex Hv1`

Check a set of gRNAs for their specificity in CRISPR/cas application. As a result for each given gRNA a score is calculated reflecting its target specififity. 

##### Use krispr with custom model

`Kmasker --run --krispr --fasta candidate_grna.fasta --kindex Hv1 --model your_model.RData`

##### Make a new model

You can create a new model for the krispr analysis. See details at [https://git.io/JecYI](https://git.io/JecYI). 

`Kmasker --run --krispr --make_model your_targets.csv`

The model will be created in the current directory. 

#### COMPARATIVE

Kmasker --run --compare --fasta query.fasta --kindex At1 Hv1

Perform a comparative study using multiple kindex. Kmasker will search for differences between both applied kindex.

#### GFF feature explanations
See also Input/Output section for a more general explanation.

Kmasker annotates so-called KRC - K-mer Repeat Clusters. A cluster consists of one or more Kmer Repeat Regions (KRR).
Sometimes regions with high repetitive contents are interrupted by regions with low repeat content or single kmers with a lower count. 
To create KRCs Kmasker first estimates seeds. Seeds can be single k-mers or larger regions with a high count.
For example, imagine the following case (X denotes a nucleotide with a high kmer count):

`XAAGTTXXXTGCXXXXXXXXXXXXXXXCCCGTAAGXXGGTXXXXXGX`

We would like to consider the whole region as a region with high kmer count. The KRC will consist of single positions or regions. 
To calculate which parts belong to the cluster and which not, the following algorithm/settings are used.  
###### Expert settings
###### mingff
Minimal length for a feature (default: KRC) to be included (and annotated) in the GFF file.
###### pctgap
The percentual gap between to regions (KRR) to be included in one cluster. Default: 10bp
The length of the gap is divided by the sum of the length of the two regions. The maximum threshold of this ratio is this parameter. Default 10% 
###### minseed
As a second calculation, there is also minseed used (if the ratio is bigger than *pctgap*). 
This is useful if we have situations like that: `XTGTTXXX`
Now the gap is just 4bp, but the sum of the regions is also just 4bp. The ratio would be 50%. But we want to include such situations in a cluster, too. 
*minseed* is the maximal distance between two (short) subfeatures (KRR) to be included in the main feature (KRC). Default: 5bp

#### Module help

```
Usage of program Kmasker:
 (version:  0.0.36 rc191017)				(session id: 72quHhdU9k)

 Command:
	 Kmasker --run --fasta sequence_to_be_analyzed.fasta


 Options:
 --kindex	 single or multiple k-mer indices (use space delimited list)
 --fasta	 sequences in FASTA format for k-mer analysis and masking
 --krispr	 provide gRNA sequences in FASTA format for specificity analysis
 --compare	 perform comparative analysis using multiple k-mer indices (requires --kindex K1 K2)
 --rept		 frequency threshold used for masking [5]!
 --minl		 minimal length of sequence. Kmasker will extract all non-repetitive sequences with sufficient length [100]
 --fish		 Extracts long sequence strechtes with low repetitiveness as FISH candidates
```


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
Link: [http://pgsb.helmholtz-muenchen.de/plant/recat/](http://pgsb.helmholtz-muenchen.de/plant/recat/)

TREP:
[http://botserv2.uzh.ch/kelldata/trep-db/](http://botserv2.uzh.ch/kelldata/trep-db/)

##### Expert settings
Kmasker blasts the KRCs against the repeat databases to estimate an annotation of the features.
You can modify the default settings of the blast commandline to better fit your needs with the expert setting --expert\_setting\_blast with the following syntax: 

`--expert_setting_blast '-perc_identity=0.95;-evalue=10'` (the default is perc_identity of 80 and evalue of 0.1)
(the syntax of the expert\_setting\_jelly also works fine: `--expert_setting_blast 'perc_identity=0.95;evalue=10'`

This example would set -perc_identity=0.95 -evalue=10 for the blast command. The parameters `-ungapped -max_hsps 1 -max_target_seqs 1` are unchangeable. The number of threads will be imported from the Kmasker `--threads` parameter.

#### Module help

```
Usage of program Kmasker:
 (version:  0.0.36 rc191017)				(session id: 71FVo0GF39)

 Command (subset):

	 Kmasker --explore --annotate --fasta query.fasta --gff kmasker_result.gff --feature KRC --dbfasta repeats.fasta

	 Kmasker --explore --hist --occ file.occ --list sequence.ids

	 Kmasker --explore --hexplot --kindex At1 Hv1

	 Kmasker --explore --stats --occ file.occ

 Options:
 --annotate		 custom annotation using featured elements of GFF (requires --gff, --fasta, --db or --dbfasta)
 --gff			 use Kmasker constructed GFF report for annotation
 --feature		 the type of feature in the GFF that should be annotated
 --dbfasta		 custom sequences [FASTA] with annotated features in sequence descriptions
 --db			 pre-calculated blastableDB of nucleotides used for annotation
 --xtract		 extract non-masked regions from Xmasked FASTA (requires --fasta)
 --hist			 create histogram using raw values (requires --occ and optional --list)
 --histm		 create histogram using calulated means (requires --occ and optional --list)
 --violin		 create violin plot (apply after '--run --compare' mode)
 --hexplot		 create hexagon plot (apply after '--run --compare' mode)
 --barplot		 create barplot (apply after '--run --compare' mode)
 --occ			 provide a Kmasker constructed occ file containing k-mer frequencies
 --file			 provide simple file (three columns) for plotting
 --cfile		 provide comparative file for plotting
 --list			 file containing a subset of selected contig identifier
 --stats		 create report of basic statistics (requires --occ and --gff)
 --cstats		 create copmarative statistics (requires --occ)
```

## INPUT/OUTPUT:

### INPUT
FASTA/FASTQ are accepted as input formats to construct KINDEX structures (Build Module). To perform k-mer analysis the user query sequence is required to be in FASTA format (Run Module). 

### OUTPUT
Various output formats are generated. Sequence output is provided in FASTA format. Detected repeats are provided in BED and GFF format. In addition, we provide statistical reports that are tab-separated files. 

#### GFF
The GFF files consist of 9 columns following standard GFF specification. Two kinds of GFF files are generated depending on the type of method that is applied. If Kmasker plants is applied for repeat detection it contains lines with KRR and KRC as feature type. If a k-mer frequency ratio analysis is performed it contains lines with KDC as feature type. KRR (‘k-mer repeat region’): Short continues nucleotide sequences with k-mer counts above threshold. KRC (‘k-mer repeat cluster’): KRR segments that are in close distance which have been merged into clusters. Merging can be adjusted with the parameter ‘--expert\_setting\_kmasker pcgap=value’. Further details are explained in the run section at expert settings.
KDC (‘k-mer diverse cluster’): Detected segments of the input sequence which show diverging k-mer patters in the comparative study of the k-mer ratio analysis. Here, the two applied sequence data sets (KINDEX A and KINDEX B) have significantly different k-mer counts.

#### OCC
The file holds the base-specific k-mer counts for one or more biological sequences stored in a corre-sponding FASTA file. It is most similar to FASTA QUAL formats. Lines starting with “>” contain the sequence identifier followed by lines with numeric values. Each values corresponds to a nucleotide position of the input sequence. K-mer counts are represented as non-negative integers separated by whitespace (typically a single space or newline), and can span multiple lines.

## Expert settings configuration files
If you have to append different expert settings for your workflow by default, you can safe them in a configuration file. 
The syntax is one parameter in each line with the value after a equal sign.
### jellyfish
append `--config_jelly=path/to/file`
e.g. `size=3`
### blast
append `--config_blast=path/to/file`
e.g.
 
```
#blast settings
perc_identity=.095
evalue=10
```
### kmasker
append `--config_kmasker=path/to/file`

The following parameters can be set in a confige file: pctgap, minseed, mingff, rept, minl and bed.


