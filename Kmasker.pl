#!/usr/bin/perl -w
use strict;
use warnings;
use IO::File;
use Getopt::Long;

my $version = "0.0.16 rc160712";
my $fasta;
my $indexfile;
my $user_index;
my $dir;
my @dirfasta;

my $build;
my $run;
my $postprocessing;
my $k_usr;
my $k = 21;
my $tool_tallymer;
my $tool_jellyfish;
my $seq_usr;
my $length_threshold	= 100;
my $length_threshold_usr;
my $repeat_threshold	= 5;
my $repeat_threshod_usr;

#GENERAL parameter
my $help;
my $keep_temporary_files;
my $show_index_repository;
my $show_list_of_species;


#DEFAULT: bowman
my $path_kindex = "/mnt/hsm_nfs/schmutzr/FREAK_KINDEXs/Barley/BOWMAN/bowman_8x/KINDEX/";
my $esa = "0dde19b3c1ac3ce3fd3f036c2e5e706e_k21";


my $result = GetOptions (	#MAIN
							"build"			=> \$build,
							"run"			=> \$run,
							"postprocessing"=> \$postprocessing,
							
							#BUILD
							"seq=s"   		=> \$seq_usr,  			# provide the fasta or fastqfile
							"tally"			=> \$tool_tallymer,
							"jelly"			=> \$tool_jellyfish,
							"k=i"			=> \$k_usr,
							
							#RUN
							"fasta=s"		=> \$fasta,							
							"pathK=s"		=> \$path_kindex,
							"rept=s"		=> \$repeat_threshod_usr,
							"min_length=s"	=> \$length_threshold_usr,
							"esa=s"			=> \$esa,							
							"repository"	=> \$show_index_repository,
							"species"		=> \$show_list_of_species,						
							"dir=s"			=> \$dir,
							"index"			=> \$user_index,
							
							#POSTPROCESSING
							"plot_hist"		=> \$plot_hist_frequency,							
							
							#GLOBAL
							"keep_tmp"		=> \$keep_temporary_files,
							"verbose"		=> \$verbose,
							"help"			=> \$help										
						);

if((defined $fasta_usr)||(defined $dir)||(defined $build)){
	# valid parameter setting 
}else{
	$help = 1;
}





if(defined $help){
	print "\n\n Usage of programm: ";
    print "\n (version:  ".$version.")";
    print "\n\n COMMAND GENERAL:";
    print "\n\n\t Kmasker --fasta sequence_to_be_analyzed.fasta\n";
    print "\n Description:\n\t Kmasker calls the whole freak pipeline automatised! \n";
    print "\t Result will contain only those region (sequences) that are longer than the \n";
    print "\t provided value (parameter 'sl') and have a lower frequency than the \n";
    print "\t threshold (parameter 'rept')";
    print "\n\n";
	print "\n --fasta\t provide FASTA input files for masking";
	print "\n --dir\t\t directory option to analyze multiple FASTA input files.";
	print "\n --build\t start construction of new index (requires --indexfiles)";
	print "\n --indexfile\t fasta sequences to build a new index";
	print "\n --index\t use specific index e.g. bowman or morex [default: bowman]";
	print "\n --pathK\t provide absolut path where all files of the index are stored";
	print "\n --esa\t\t core name of index";
	print "\n --rept\t\t frequency threshold used for masking [5]!";
	print "\n --sl\t\t sequence length threshold. After masking sequences have to be longer than provided length threshold.";
	print "\n\n";
	print "\n COMMAND RUN pre-caculated barley index [bowman or morex] :";
	print "\n";
	print "\n\t >Kmasker --fasta sequence_to_be_analyzed.fasta --index morex";
	print "\n";
	print "\n COMMAND BUILD index :";
	print "\n";
	print "\n\t >Kmasker --build --seq mysequences.fasta";
	print "\n";	
	print "\n COMMAND APPLY index (using local custom index):";
	print "\n";
	print "\n\t >Kmasker --fasta sequence_to_be_analyzed.fasta --pathK /abolute/path/to/index --esa index.name";
	print "\n\n";
	exit();
}

#details on the index and the default frequency level:
# Kmasker is working with the parameter rept that sets the real frequency. 
# Since the KINDEX used here is based on 8x Barley sequencing, the real 
# frequency has to be devided by eight! Thus if you use e.g. rept 40 in freak it
# is a real frequency in the genome of 5!
#


1;
