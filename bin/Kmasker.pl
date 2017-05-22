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
use kmasker::kmasker_build qw(build_kindex_jelly make_config remove_repository_entry set_kindex_global clean_repository_directory set_private_path);
use kmasker::kmasker_run qw(run_kmasker_SK run_kmasker_MK show_version_PM_run);
use kmasker::kmasker_postprocessing qw(plot_histogram);

my $version 	= "0.0.24 rc170506";
my $path 		= dirname abs_path $0;		
my $fasta;
my $fastq;
my $indexfile;

#MODULES
#BUILD
my $build;
my $run;
my $postprocessing;
my $repositories;
my $build_config;
my $make_config;
my $genome_size = "N/A";
my $genome_size_usr;
my $name;
my $name_usr;
my $PATH_kindex_private = "";
my $PATH_kindex_global 	= "";

#RUN
my $kindex_usr;
my $k_usr;
my $k 					= 21;
my $tool_jellyfish;
my @seq_usr;
my $length_threshold	= 100;
my $length_threshold_usr;
my $repeat_threshold	= 5;
my $repeat_threshod_usr;
my $tolerant_length_threshold_usr;
my $tolerant_length_threshold = 0;
my $MK_percent_gapsize	= 10;	#default	FIXME
my $MK_min_seed			= 5;	#default	FIXME
my $MK_min_gff			= 10;   #default	FIXME

#Postprocessing
my $gff;
my $repeat_lib_user;
my $repeat_lib			= "REdat";
my $clist;
my $occ;
my $stats;

#GENERAL parameter
my $help;
my $keep_temporary_files;
my $show_kindex_repository;
my $show_details_for_kindex;
my $set_private_path;
my $remove_kindex;
my $plot_hist_frequency;
my $expert_setting = ""; 
my $set_global;
my $user_name;
my $verbose;

#HASH
my %HASH_repository_kindex;
my %HASH_path;

#DEFAULT: no default anymore
my $kindex;
my @multi_kindex;

my $result = GetOptions (	#MAIN
							"build"				=> \$build,
							"run"				=> \$run,
							"postprocessing"	=> \$postprocessing,
							
							#BUILD
							"seq=s{1,}"   		=> \@seq_usr,  			# provide the fasta or fastqfile
							"k=i"				=> \$k_usr,
							"gs=i"				=> \$genome_size_usr,
							"config=s"			=> \$build_config,
							"make_config"		=> \$make_config,
							"name=s"			=> \$name_usr,
							
							#RUN
							"fasta=s"			=> \$fasta,	
							"kindex=s"			=> \$kindex_usr,
							"multi_kindex=s{1,}"=> \@multi_kindex,
							"rept=s"			=> \$repeat_threshod_usr,
							"min_length=s"		=> \$length_threshold_usr,
#							"tol_length=s"		=> \$tolerant_length_threshold_usr,
												
							#POSTPROCESSING
							"plot_hist"			=> \$plot_hist_frequency,
							"clist=s"			=> \$clist,
							"occ=s"				=> \$occ,
							"stats"				=> \$stats,
#							"gff"				=> \$gff,
#							"repeat_library=s"	=> \$repeat_lib_user,							
							
							#GLOBAL
							"show_repository"	=> \$show_kindex_repository,
							"show_details=s"	=> \$show_details_for_kindex,
							"remove_kindex=s"	=> \$remove_kindex,
							"set_global=s"		=> \$set_global,
							"set_private_path=s"=> \$set_private_path,
							
							#Houskeeping
							"expert_setting"	=> \$expert_setting,
							"keep_tmp"			=> \$keep_temporary_files,
							"verbose"			=> \$verbose,
							"help"				=> \$help										
						);
						


if(defined $help){
	
	print "\n Usage of program Kmasker: ";
    print "\n (version:  ".$version.")";
    print "\n";
	
	if(defined $build){
		#HELP section build		
		print "\n Command:";
		print "\n\t Kmasker --build --seq mysequences.fasta";
		
		print "\n\n Options:";
		print "\n --seq\t\t fasta or fastq sequence(s) that are used to build the index";
		print "\n --k\t\t k-mer size to build index [21]";
		print "\n --gs\t\t genome size of species (in Mbp)";
		print "\n --name \t\t provide name to identify the kindex (e.g. HvMRX for hordeum vulgare cultivare morex) [date]";
		print "\n --make_config\t creates basic config file ('build_kindex.config') for completion by user";
		print "\n --config\t configuration file providing information used for construction of kindex";
		
		print "\n\n";
		exit();
	}
	
	if(defined $run){
		#HELP section run
		print "\n Command:";
		print "\n\t Kmasker --run --fasta sequence_to_be_analyzed.fasta\n";		
		
		print "\n\n Options:";
		print "\n --fasta\t FASTA sequence for k-mer analysis and masking";
		print "\n --kindex\t use specific k-mer index e.g. bowman or morex";
		print "\n --multi_kindex\t use multiple k-mer indices for comparative analysis of FASTA sequence (e.g. bowman and morex)";
		print "\n --rept\t\t frequency threshold used for masking [5]!";
		print "\n --min_length\t minimal length of sequence. Kmasker will extract all non-repetitive sequences with sufficient length [100]";
#		print "\n --tol_length\t maximal length of sequence with high k-mer frequencies. Within non-repetitive candidate sequences with sufficient sequence length \
#		 		\t\t\t	(--min_length) it is tolerated that small regions occure were the corresponding k-mer frequency exceeds the defined threshold (--rept). [0]";
	
		print "\n\n";
		exit();
	}
	
	if(defined $postprocessing){
		#HELP section postprocessing
		print "\n Command:";
		print "\n\t Kmasker --postprocessing --plot_history --occ file.occ --clist list_of_contigs.txt";
		
		print "\n\n Options:";
		print "\n --occ\t\t provide a Kmasker constructed occ file containing k-mer frequencies";
		print "\n --plot_hist\t\t create graphical output as histogram (requires --clist)";
		print "\n --clist\t\t file containing a list of contig identifier that are used in postprocessing";	

#		print "\n --stats\t\t\t calculate basic statistics like avegare k-mer frequency per contig etc. (requires --occ)";	
#		print "\n --gff\t\t\t perform repeat annotation and construct GFF report";
#		print "\n --repeat_library\t provide repeat library [REdat]"; 
		
		print "\n\n";
		exit();
	}
	

    print "\n Description:\n\t Kmasker is a tool for the automatic detection of repetitive sequence regions.";
    
    print "\n\n Modules:";
	print "\n --build\t\t construction of new index (requires --indexfiles)";
	print "\n --run\t\t\t run k-mer repeat detection and masking (requires --fasta)";
	print "\n --postprocessing\t perform downstream analysis with constructed index and detected repeats";
	
	print "\n\n General options:";
	print "\n --show_repository\t shows complete list of global and private k-mer indices";
	print "\n --show_details\t\t shows details for a requested kindex";
	print "\n --remove_kindex\t remove kindex from repository";
	print "\n --expert_setting\t submit individual parameter to Kmasker (e.g. on memory usage for index construction)";
	
	print "\n\n";
	exit();
}


##MAIN

#load global settings
&read_user_config;
&read_user_repositories;

#USER specification
#kindex
if(defined $kindex_usr){
	if(exists $HASH_repository_kindex{$kindex_usr}){
		$kindex = $kindex_usr;
	}else{
		print "\n ERROR: defined kindex ('".$kindex_usr."') does not exist!\n\n";
		exit();
	}	
}	

#name
if(defined $name_usr){
	$name = $name_usr;
}

#genome_size
if(defined $genome_size_usr){
	if($genome_size_usr =~ /^[+-]?\d+$/){
		#is number
		$genome_size = $genome_size_usr;
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

#CHECK setting
&check_settings;

if(defined $build){
	#USE BUILD MODULE
	
	if(defined $make_config){
		&make_config;
	}elsif(defined $set_private_path){
		&set_private_path($set_private_path, \%HASH_repository_kindex);
	}else{
		#INIT
		my %HASH_info 						= ();
		my $input 							= join(" ", sort { $a cmp $b } @seq_usr);
		$HASH_info{"user_name"}				= $user_name;
		$HASH_info{"seq"} 					= $input;
		$HASH_info{"k-mer"}					= $k;
		$HASH_info{"genome_size"}			= $genome_size;
		$HASH_info{"short_tag"}				= $name;
		$HASH_info{"expert_setting"}		= $expert_setting;
		$HASH_info{"PATH_kindex_global"}	= $PATH_kindex_global; 
		$HASH_info{"PATH_kindex_private"}	= $PATH_kindex_private;
		$HASH_info{"path_bin"}				= $path; 
		
		#EDIT
		if(defined $set_global){
			if(exists $HASH_repository_kindex{$set_global}){
				&set_kindex_global($set_global, \%HASH_info);			
			}else{
				print "\n WARNING: requested kindex '".$set_global."' does not exist. Kmasker was stopped!\n\n";
			}
			exit();			
		}
		
		#CONSTRUCT
		if(defined $input){
			&build_kindex_jelly(\%HASH_info, $build_config, \%HASH_repository_kindex, \%HASH_path); 
		}
		
		#CLEAN
		&read_user_repositories;
		&clean_repository_directory(\%HASH_info, \%HASH_repository_kindex);			
	}	
		
	#QUIT
	print "\n - Thanks for using Kmasker! -\n\n";
	exit();
}

if(defined $run){
	#USE RUN MODULE
	
	my %HASH_info 						= ();
	$HASH_info{"user_name"}				= $user_name;
	$HASH_info{"rept"}					= $repeat_threshold; 
	$HASH_info{"MK_percent_gapsize"}	= $MK_percent_gapsize;
	$HASH_info{"MK_min_seed"}			= $MK_min_seed;
	$HASH_info{"MK_min_gff"}			= $MK_min_gff;
	$HASH_info{"expert_setting"}		= $expert_setting; 
	$HASH_info{"PATH_kindex_global"}	= $PATH_kindex_global; 
	$HASH_info{"PATH_kindex_private"}	= $PATH_kindex_private;
	$HASH_info{"path_bin"}				= $path;
	
	if(defined $kindex){
	#single kindex				
		&run_kmasker_SK($fasta, $kindex, \%HASH_info, \%HASH_repository_kindex);
	}elsif(scalar(@multi_kindex > 1)){
	#multiple kindex
		&run_kmasker_MK($fasta, \@multi_kindex, \%HASH_info, \%HASH_repository_kindex);
	}
	
	#QUIT
	print "\n - Thanks for using Kmasker! -\n\n";
	exit();
}

if(defined $postprocessing){
	#USE POSTPROCESSING MODULE
	
	if(defined $occ){
	# postprocessing requires an OCC file
		
		my $missing_parameter = "";
	
		if(defined $plot_hist_frequency){
			if(defined $clist){
				&plot_histogram($occ, $clist);
			}else{
				$missing_parameter .= " --clist";
			}
		}		
		
		if($missing_parameter ne ""){
			#GIVE warning note for missing parameter
			print "\n ERROR: missing parameter (".$missing_parameter.") !\n\n";
		}
		
	}else{
		print "\n ERROR: no occ provided. For Kmasker postprocessing an occ file is required!\n\n";
	}
	
	#QUIT
	print "\n - Thanks for using Kmasker! -\n\n";
	exit();
}
	

#GENERAL options
if(defined $show_kindex_repository){
	&show_repository();
	exit();
}	
	
if(defined $show_details_for_kindex){
	&show_details_for_kindex($show_details_for_kindex);
	exit();
}

if(defined $remove_kindex){
	my %HASH_info 						= ();
	$HASH_info{"user_name"}				= $user_name;
	$HASH_info{"PATH_kindex_global"}	= $PATH_kindex_global; 
	$HASH_info{"PATH_kindex_private"}	= $PATH_kindex_private; 
	&remove_repository_entry($remove_kindex,\%HASH_info);
	exit();
}				

##END MAIN


## subroutine
#
sub read_user_config(){
	$user_name 			= `whoami`;
	$user_name			=~ s/\n//g;
	my $gconf 			= $path."/config.kmasker";
	my $uconf 			= $ENV{"HOME"}."/.user_config.kmasker";
	my $urepositories 	= $ENV{"HOME"}."/.user_repositories.kmasker";
	
	if(-e $gconf){
		#LOAD global info
		my $gCFG = new IO::File($gconf, "r") or die "\n unable to read user config $!";	
		
		while(<$gCFG>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_tmp = split("\t", $line);
			$PATH_kindex_global	= $ARRAY_tmp[1] if($ARRAY_tmp[0] eq "PATH_kindex_global");
			$PATH_kindex_global .= "/" if($PATH_kindex_global !~ /\/$/);
			
			#GLOBAL
			#PRIVATE
			if(-d $PATH_kindex_global){
				#directory exists - do nothing
			}else{
				#directory has to be created
				system("mkdir ".$PATH_kindex_global);
			}
			
			#READ external tool path
			#JELLYFISH
			if($line =~ /^jellyfish=/){
				my @ARRAY_tmp = split("=", $line);
				if(!defined $ARRAY_tmp[1]){
					system("command -v jellyfish >/dev/null 2>&1 || { echo >&2 \"Kmasker requires jellyfish but it's not installed! Kmasker process stopped.\"; exit 1; \}");
					$HASH_path{"jellyfish"} = `command -v jellyfish`;
					$HASH_path{"jellyfish"} =~ s/\n//;
				}else{
					$HASH_path{"jellyfish"} = $ARRAY_tmp[1];
					system("command -v ".$HASH_path{"jellyfish"}." >/dev/null 2>&1 || { echo >&2 \"Kmasker requires jellyfish but it's not installed!  Kmasker process stopped.\"; exit 1; \}");
				}
				print "\n jellyfish=".$HASH_path{"jellyfish"}."\n" if(defined $verbose);
			}
			
			#FASTQ-STATs
			if($line =~ /^fastq-stats=/){
				my @ARRAY_tmp = split("=", $line);
				if(!defined $ARRAY_tmp[1]){
					system("command -v fastq-stats >/dev/null 2>&1 || { echo >&2 \"Kmasker requires fastq-stats but it's not installed! Kmasker process stopped.\"; exit 1; \}");
					$HASH_path{"fastq-stats"} = `command -v fastq-stats`;
					$HASH_path{"fastq-stats"} =~ s/\n//;
				}else{
					$HASH_path{"fastq-stats"} = $ARRAY_tmp[1];
					system("command -v ".$HASH_path{"fastq-stats"}." >/dev/null 2>&1 || { echo >&2 \"Kmasker requires fastq-stats but it's not installed! Kmasker process stopped.\"; exit 1; \}");
				}
				print "\n fastq-stats=".$HASH_path{"fastq-stats"}."\n" if(defined $verbose);
			}
			
			#GFFREAD
			if($line =~ /^gffread=/){
				my @ARRAY_tmp = split("=", $line);
				if(!defined $ARRAY_tmp[1]){
					system("command -v gffread >/dev/null 2>&1 || { echo >&2 \"Kmasker requires gffread but it's not installed! Kmasker process stopped.\"; exit 1; \}");
					$HASH_path{"gffread"} = `command -v gffread`;
					$HASH_path{"gffread"} =~ s/\n//;
				}else{
					$HASH_path{"gffread"} = $ARRAY_tmp[1];
					system("command -v ".$HASH_path{"gffread"}." >/dev/null 2>&1 || { echo >&2 \"Kmasker requires gffread but it's not installed! Kmasker process stopped.\"; exit 1; \}");
				}
				print "\n gffread=".$HASH_path{"gffread"}."\n" if(defined $verbose);
			}
		}
	}
	
	if(-e $uconf){
		#LOAD private info
		my $uCFG = new IO::File($uconf, "r") or die "\n unable to read user config $!";	
		while(<$uCFG>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_tmp = split("\t", $line);
			$PATH_kindex_private= $ARRAY_tmp[1] if($ARRAY_tmp[0] eq "PATH_kindex_private");
			$PATH_kindex_private.= "/" if($PATH_kindex_private !~ /\/$/);
			
			#PRIVATE
			if(-d $PATH_kindex_private){
				#directory exists - do nothing
			}else{
				#directory has to be created
				system("mkdir ".$PATH_kindex_private);
			}			
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
	my $urepositories 	= $ENV{"HOME"}."/.user_repositories.kmasker";
	if(-e $urepositories){
		#LOAD info
		
		#load global repository
		my $REPO_global = new IO::File($path."/repositories.kmasker", "r") or die "\n unable to read global repository list $!";	
		while(<$REPO_global>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_line = split("\t", $line);
			$HASH_repository_kindex{$ARRAY_line[0]} = $line;
		}
		
		#load priavte user repositories
		my $REPO_private = new IO::File($ENV{"HOME"}."/.user_repositories.kmasker", "r") or die "\n unable to read user repository list $!";	
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
sub show_repository(){
	
	print "\n\nREPOSITORY of available kindex structures:\n";
	
	foreach my $kindex_this (keys %HASH_repository_kindex){
		my @ARRAY_line = split("\t", $HASH_repository_kindex{$kindex_this});
		print "\n\t".$ARRAY_line[0]."\t\t".$ARRAY_line[1];
		if(defined $ARRAY_line[11]){
			print "\t".$ARRAY_line[11];
		}
	}
	print "\n\n";
}


## subroutine
#
sub show_details_for_kindex(){
	my $kindex = $_[0];
	if(exists $HASH_repository_kindex{$kindex}){
		my $linearray = $HASH_repository_kindex{$kindex};
		my @ARRAY_details = split("\t", $linearray);
		
		print "\n\n KINDEX details for ".$kindex." \n";
		print "\n\tcommon_name:      ".$ARRAY_details[1];
		print "\n\tscientific_name:  ".$ARRAY_details[2];
		if(defined $ARRAY_details[12]){
			print "\n\tgenome_size       ".$ARRAY_details[12]." [Mbp]";
		}else{
			print "\n\tgenome_size       N/A";
		}
		print "\n\ttype              ".$ARRAY_details[3];
		print "\n\tsequencing_depth: ".$ARRAY_details[4];
		print "\n\tk-mer:            ".$ARRAY_details[5];
		print "\n\tsequence_type:    ".$ARRAY_details[6];
		print "\n\tnote:             ".$ARRAY_details[7];
		print "\n\tmd5sum:           ".$ARRAY_details[8];
		print "\n\tseq:              ".$ARRAY_details[9];
		print "\n\tabsolut_path:     ".$ARRAY_details[10];
		print "\n\tstatus:           ".$ARRAY_details[11];	
		print "\n\n";
	}else{
		print "\n\n WARNING: Requested kindex (".$kindex."). does not exist. Please check and use different index name.\n\n";
		exit();
	}
}


## subroutine
#
sub initiate_user(){
	
	$user_name 			= `whoami`;
	$user_name			=~ s/\n//g;
	my $uconf 			= $ENV{"HOME"}."/.user_config.kmasker";
	my $urepositories 	= $ENV{"HOME"}."/.user_repositories.kmasker";
	if(-e $uconf){
		#USER already exists, do nothing
	}else{
		#SETUP user conf
		my $USER_CONF 	= new IO::File($uconf, "w") or die "could not write user repository : $!\n";
		print $USER_CONF "PATH_kindex_private\t".$ENV{"HOME"}."/KINDEX/";
		close $USER_CONF;	
		print "\n PLEASE NOTE: \n You are writing all large data structures to your home directory [default].";
		print "\n It is recommended to modify the path for 'PATH_kindex_private'.\n";
		print "\n Use the following command: 'Kmasker --build --set_private_path enter/your/path'\n\n";	
		
		#SETUP user repo
		my $USER_REPOS 	= new IO::File($urepositories, "w") or die "could not write user repository : $!\n";
		
		print $USER_REPOS "##KINDEX repository\n";
		#1
		print $USER_REPOS "#short_tag\t";
		#2
		print $USER_REPOS "common_name\t";
		#3
		print $USER_REPOS "scientific_name\t";
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
		#11
		print $USER_REPOS "status\n";	
	}	
}

## subroutine
#
sub check_settings(){
	
	if($PATH_kindex_private eq $ENV{"HOME"}."/KINDEX/"){
		print "\n PLEASE NOTE: \n You are writing all large data structures to your home directory [default].";
		print "\n It is recommended to modify the path for 'PATH_kindex_private'.\n";
		print "\n Use the following command: 'Kmasker --build --set_private_path enter/your/path'\n\n";
	}
	
	my $module_count = 0;
	if(defined $build){
		$module_count++;
		
		if(defined $make_config){
			#nothing to do
		}elsif(defined $set_global){
			#nothing to do
		}elsif(defined $set_private_path){
			#nothing to do
		}elsif((scalar @seq_usr) == 0){
			print "\n .. kmasker was stopped: no input sequence provided (--seq) !";
			print "\n\n";
			exit(0);
		}
	}
	if(defined $run){
		$module_count++;
		if(!((defined $kindex)||(scalar (@multi_kindex >1)))){
			print "\n .. kmasker was stopped: no kindex defined (--kindex) !";
			print "\n\n";
			exit(0);
		}
		
		if(!(defined $fasta)){
			print "\n .. kmasker was stopped: no sequence provided (--fasta) !";
			print "\n\n";
			exit(0);
		}
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
