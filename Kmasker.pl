#!/usr/bin/perl -w
use strict;
use warnings;
use IO::File;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex md5_base64);

#setup package directory
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname abs_path $0) . '/lib';

#include packages
use kmasker::kmasker_build qw(build_kindex_jelly make_config);

my $version = "0.0.17 rc160929";
#my $CORE 	= "/opt/Bio/Kmasker/0.0.15/bin";
my $CORE 	= "/data/filer/agbi/schmutzr/projects/KMASKER/source_code_v160929_Kmasker_plants";

my $fasta;
my $fastq;
my $indexfile;
my $user_kindex;

#MODULES
my $build;
my $run;
my $postprocessing;
my $repositories;
my $build_config;
my $make_config;
my $PATH_kindex_private = "";
my $PATH_kindex_global 	= "";

my $kindex_usr;
my $k_usr;
my $k 					= 21;
my $tool_jellyfish;
my @seq_usr;
my $length_threshold	= 100;
my $length_threshold_usr;
my $repeat_threshold	= 5;
my $repeat_threshod_usr;

#Postprocessing
my $gff;
my $repeat_lib_user;
my $repeat_lib			= "REdat";

#GENERAL parameter
my $help;
my $keep_temporary_files;
my $show_index_repository;
my $show_list_of_species;
my $plot_hist_frequency;
my $user_name;
my $verbose;

#HASH
my %HASH_repository_kindex;
my %HASH_repository_kindex_short_tag;

#DEFAULT: bowman
my $kindex 		= "Hv8x";

my $result = GetOptions (	#MAIN
							"build"			=> \$build,
							"run"			=> \$run,
							"postprocessing"=> \$postprocessing,
							
							#BUILD
							"seq=s{1,}"   	=> \@seq_usr,  			# provide the fasta or fastqfile
							"k=i"			=> \$k_usr,
							"config=s"		=> \$build_config,
							"make_config"	=> \$make_config,
							
							#RUN
							"fasta=s"		=> \$fasta,	
							"kindex=s"		=> \$kindex_usr,
							"rept=s"		=> \$repeat_threshod_usr,
							"min_length=s"	=> \$length_threshold_usr,
							"index"			=> \$user_kindex,
							"repositories"	=> \$show_index_repository,
							"species"		=> \$show_list_of_species,
							
							#POSTPROCESSING
							"plot_hist"		=> \$plot_hist_frequency,	
							"gff"			=> \$gff,
							"repeat_library=s"	=> \$repeat_lib_user,							
							
							#GLOBAL
							"keep_tmp"		=> \$keep_temporary_files,
							"verbose"		=> \$verbose,
							"help"			=> \$help										
						);
						


if(defined $help){
	
	print "\n Usage of program Kmasker: ";
    print "\n (version:  ".$version.")";
    print "\n";
	
	if(defined $build){
		#HELP section build		
		print "\n Command:";
		print "\n\t Kmasker --build --seq mysequences.fasta";
		
		print "\n\n Option(s):";
		print "\n --seq\t\t fasta or fastq sequence(s) that are used to build the index";
		print "\n --k\t\t k-mer size to build index [21]";
		print "\n --make_config\t creates basic config file ('build_kindex.config') for completion by user";
		print "\n --config\t configuration file providing information used for construction of kindex";
		
		print "\n\n";
		exit();
	}
	
	if(defined $run){
		#HELP section run
		print "\n Command:";
		print "\n\t Kmasker --fasta sequence_to_be_analyzed.fasta\n";		
		
		print "\n\n Option(s):";
		print "\n --fasta\t FASTA sequence for k-mer analysis and masking";
		print "\n --kindex\t use specific k-mer index e.g. bowman or morex [default: bowman]";
		print "\n --multi_kindex\t use multiple k-mer indices for comparative analysis of FASTA sequence (e.g. bowman and morex)";
		print "\n --rept\t\t frequency threshold used for masking [5]!";
		print "\n --min_length\t minimal length of sequence. Kmasker will extract all non-repetitive sequences with sufficient length [100]";
		print "\n --repositories\t show complete list of global and private repositories of existing k-mer indices";
		print "\n --species\t show list of existing species with build kindex";		
	
		print "\n\n";
		exit();
	}
	
	if(defined $postprocessing){
		#HELP section postprocessing
		print "\n Command:";
		print "\n\t Kmasker --postprocessing --gff mysequences.fasta";
		
		print "\n\n Option(s):";
		print "\n --plot_hist\t\t create graphical output as histogram (looks for --clist)";
		print "\n --clist\t\t list of contig identifier that are used in postprocessing";		
		print "\n --gff\t\t\t perform repeat annotation and construct GFF report";
		print "\n --repeat_library\t provide repeat library [REdat]"; 
		
		print "\n\n";
		exit();
	}
	

    
    print "\n Description:\n\t Kmasker is a tool for the automatic detection of repetitive sequence regions.";
    
    print "\n\n Option(s):";
	print "\n --build\t\t construction of new index (requires --indexfiles)";
	print "\n --run\t\t\t run k-mer repeat detection and masking (requires --fasta)";
	print "\n --postprocessing\t perform downstream analysis with constructed index and detected repeats";
	print "\n\n";
	exit();
}


##MAIN

#load global settings
&check_settings;
&read_user_config;
&read_user_repositories;

#USER specification
#kindex
if(defined $user_kindex){
	if(exists $HASH_repository_kindex_short_tag{$user_kindex}){
		$kindex = $user_kindex;
	}else{
		print "\n ERROR: defined kindex ('".$user_kindex."') does not exist!";
		exit();
	}	
}	

#k-mer size
if(defined $k_usr){
	if($k_usr =~ /^[+-]?\d+$/){
		#is number
		$k = $k_usr;
	}
}	

#rept
if(defined $repeat_threshod_usr){
	if($repeat_threshod_usr =~ /^[+-]?\d+$/){
		#is number
		$repeat_threshold = $repeat_threshod_usr;
	}
}

#min length
if(defined $length_threshold_usr){
	if($length_threshold_usr =~ /^[+-]?\d+$/){
		#is number
		$length_threshold = $length_threshold_usr;
	}
}

#repeat library
if(defined $repeat_lib_user	){
	#FIX
}

if(defined $build){
	#USE BUILD MODULE
	
	if(defined $make_config){
		&make_config;
	}else{
	
		my %HASH_info 						= ();
		my $input 							= join(" ", sort { $a cmp $b } @seq_usr);
		$HASH_info{"user_name"}				= $user_name;
		$HASH_info{"seq"} 					= $input;
		$HASH_info{"k-mer"}					= $k; 
		$HASH_info{"PATH_kindex_global"}	= $PATH_kindex_global; 
		$HASH_info{"PATH_kindex_private"}	= $PATH_kindex_private; 
		
		&build_kindex_jelly(\%HASH_info, $build_config, \%HASH_repository_kindex); 
	}
}

if(defined $run){
	#USE BUILD MODULE
	&run_kmasker($fasta, $kindex, \%HASH_repository_kindex);	 
}

if(defined $postprocessing){
	#USE BUILD MODULE
	 
}				

##END MAIN


## subroutine
#
sub read_user_config(){
	$user_name 			= `whoami`;
	$user_name			=~ s/\n//g;
	my $uconf 			= "/home/".$user_name."/.user_config.kmasker";
	my $urepositories 	= "/home/".$user_name."/.user_repositories.kmasker";
	if(-e $uconf){
		#LOAD info
		my $uCFG = new IO::File($uconf, "r") or die "\n unable to read user config $!";	
		while(<$uCFG>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_tmp = split("\t", $line);
			$PATH_kindex_global	= $ARRAY_tmp[1] if($ARRAY_tmp[0] eq "PATH_kindex_global");
			$PATH_kindex_global .= "/" if($PATH_kindex_global !~ /\/$/);			
			$PATH_kindex_private= $ARRAY_tmp[1] if($ARRAY_tmp[0] eq "PATH_kindex_private");
			$PATH_kindex_private.= "/" if($PATH_kindex_private !~ /\/$/);
		}
	}else{
		#SETUP user
		&initiate_user();
	}
}

## subroutine
#
sub read_user_repositories(){
	#if /home/user/.user_config.kmasker not exist --> initiate
	
	$user_name 			= `whoami`;
	$user_name			=~ s/\n//g;
	my $urepositories 	= "/home/".$user_name."/.user_repositories.kmasker";
	if(-e $urepositories){
		#LOAD info
		
		#load global repository
		my $REPO_global = new IO::File($CORE."/bin/repositories.kmasker", "r") or die "\n unable to read global repository list $!";	
		while(<$REPO_global>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_line = split("\t", $line);
			$HASH_repository_kindex{$ARRAY_line[0]} = $line;
		}
		
		#load priavte user repositories
		my $REPO_private = new IO::File("/home/".$user_name."/.user_repositories.kmasker", "r") or die "\n unable to read user repository list $!";	
		while(<$REPO_private>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_line = split("\t", $line);
			$HASH_repository_kindex{$ARRAY_line[0]} = $line;
		}
		
	}else{
		#SETUP user
		&initiate_user();
	}	
}

## subroutine
#
sub initiate_user(){
	
	$user_name 			= `whoami`;
	$user_name			=~ s/\n//g;
	my $uconf 			= "/home/".$user_name."/.user_config.kmasker";
	my $urepositories 	= "/home/".$user_name."/.user_repositories.kmasker";
	if(-e $uconf){
		#USER already exists, do nothing
	}else{
		#SETUP user
		system("cp ".$CORE."/bin/config.kmasker /home/".$user_name."/.user_config.kmasker");
		my $USER_REPOS 	= new IO::File("/home/".$user_name."/.user_repositories.kmasker", "w") or die "could not write user repository : $!\n";
		
		print $USER_REPOS "##KINDEX repository\n";
		#1
		print $USER_REPOS "#short_tag\t";
		#2
		print $USER_REPOS "species\t";
		#3
		print $USER_REPOS "botanic_name\t";
		#4
		print $USER_REPOS "type (cultivar;genotype;etc.)\t";
		#5
		print $USER_REPOS "sequencing_depth\t";
		#6
		print $USER_REPOS "k-mer\t";
		#7
		print $USER_REPOS "sequence_type\t";
		#8
		print $USER_REPOS "note\t";
		#9
		print $USER_REPOS "md5sum (filename)\t";
		#10
		print $USER_REPOS "absolut_path (filename)\n";		
	}	
}

## subroutine
#
sub check_settings(){
	
	my $module_count = 0;
	if(defined $build){
		$module_count++;
		
		if(defined $make_config){
			#nothing to do
		}elsif((scalar @seq_usr) == 0){
			print "\n .. kmasker was stopped: no input sequence provided (--seq) !";
			print "\n\n";
			exit(0);
		}
	}
	if(defined $run){
		$module_count++;
	}
	if(defined $postprocessing){
		$module_count++; 
	}
	
	#MULTIPLE 
	if($module_count > 1){
		print "\n Kmasker was stopped. Multiple modules (build, run or postprocessing were used!\n";
		exit(0);
	}	
}




1;
