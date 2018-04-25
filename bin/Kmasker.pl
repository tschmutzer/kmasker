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
use kmasker::kmasker_build qw(build_kindex_jelly remove_kindex set_kindex_global set_private_path set_external_path clean_repository_directory read_config);
use kmasker::kmasker_run qw(run_kmasker_SK run_kmasker_MK show_version_PM_run);
use kmasker::kmasker_explore qw(plot_histogram custom_annotation);

my $version 	= "0.0.29 rc180425";
my $path 		= dirname abs_path $0;		
my $fasta;
my $indexfile;

#MODULES
#BUILD
my $build;
my $run;
my $explore;
my $repositories;
my $build_config;
my $make_config;
my $genome_size;
my $genome_size_usr;
my $common_name 		= "";
my $common_name_usr ;
my $index_name;
my $index_name_usr;
my $threads_usr;
my $threads = 4;
my $mem_usr;
my $mem		= 24;
my $PATH_kindex_private = "";
my $PATH_kindex_global 	= "";
my $PATH_kindex_external= "";

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
# set the following parameter using expert_setting_kmasker
my $MK_percent_gapsize	= 10;	#Default	MK_percent_gapsize (N) is the paramater allowing N percent of a region to be not repeatitive when merging regions. Used in SK and MK.
my $MK_min_seed			= 5;	#Default	MK_min_seed (N) is the minimal length (bp) of a featuere to start in the process of merging regions for GFF
my $MK_min_gff			= 10;   #Default	MK_min_gff (N) is the minimal length (bp) of a featue to be reported in GFF

#EXPLORE
my $gff;
my $clist;
my $occ;
my $stats;
my $custom_annotate;
my $blastableDB;
my $dbfasta;
#visualisation
my $hist;
my $violin;
my $hexplot;
my $triplot;

#GENERAL parameter
my $help;
my $keep_temporary_files;
my $show_kindex_repository;
my $show_details_for_kindex;
my $set_private_path;
my $set_global_path;
my $set_external_path;
my $check_install;
my $remove_kindex;
my $expert_setting_kmasker 	= "",
my $expert_setting_jelly 	= ""; 
my $expert_setting_blast	= "";
my $set_global;
my $user_name;
my $verbose;
my $temp_path			= "./temp/";

#HASH
my %HASH_repository_kindex;
my %HASH_path;

#DEFAULT: no default anymore
my $kindex;
my @multi_kindex;

my $result = GetOptions (	#MAIN
							"build"				=> \$build,
							"run"				=> \$run,
							"explore"			=> \$explore,
							
							#BUILD
							"seq=s{1,}"   		=> \@seq_usr,  			# provide the fasta or fastqfile
							"k=i"				=> \$k_usr,
							"gs=i"				=> \$genome_size_usr,
							"cn=s"				=> \$common_name_usr,
							"in=s"				=> \$index_name_usr,
							"config=s"			=> \$build_config,
							"make_config"		=> \$make_config,
							
							
							#RUN
							"fasta=s"			=> \$fasta,	
							"kindex=s"			=> \$kindex_usr,
							"multi_kindex=s{1,}"=> \@multi_kindex,
							"rept=s"			=> \$repeat_threshod_usr,
							"min_length=s"		=> \$length_threshold_usr,
#							"tol_length=s"		=> \$tolerant_length_threshold_usr,

								
												
							#EXPLORE
							"annotate"			=> \$custom_annotate,
							"gff=s"				=> \$gff,
							"db=s"				=> \$blastableDB,
							"dbfasta=s"			=> \$dbfasta,		
							"hexplot"			=> \$hexplot,
							"hist"				=> \$hist,
							"violin"			=> \$violin,
							"triplot"			=> \$triplot,
							"clist=s"			=> \$clist,
							"occ=s"				=> \$occ,
							"stats"				=> \$stats,
							
							
							#GLOBAL
							"show_repository"		=> \$show_kindex_repository,
							"show_details=s"		=> \$show_details_for_kindex,
							"remove_kindex=s"		=> \$remove_kindex,
							"set_global=s"			=> \$set_global,
							"set_private_path=s"	=> \$set_private_path,
							"set_external_path=s"	=> \$set_external_path,
							"check_install"			=> \$check_install,	
							"threads=i"				=> \$threads_usr,
							"mem=i"					=> \$mem_usr,
							
							#Houskeeping
							"expert_setting_jelly"	=> \$expert_setting_jelly,
							"expert_setting_kmasker"=> \$expert_setting_kmasker,
							"keep_tmp"				=> \$keep_temporary_files,
							"verbose"				=> \$verbose,
							"help"					=> \$help										
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
		print "\n --in \t\t provide k-mer index name (e.g. HvMRX for hordeum vulgare cultivare morex) [date]";
		print "\n --cn \t\t provide common name of species (e.g. barley)";
		print "\n --config\t configuration file providing information used for construction of kindex";
		print "\n --mem\t\t set memory limit in gigabyte [24]";
		
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
	
	if(defined $explore){
		#HELP section explore
		print "\n Command (subset):";
		print "\n\n\t Kmasker --explore --annotate --fasta query.fasta --gff kmasker_result.gff --dbfasta repeats.fasta";
		print "\n\n\t Kmasker --explore --hist --occ file.occ --clist list_of_contigs.txt";
		print "\n\n\t Kmasker --explore --hexplot --multi_kindex At1 Hv1";
		
		print "\n\n Options:";
		print "\n --annotate\t\t custom annotation using featured elements of GFF (requires --gff, --fasta, --db or --dbfasta)";
		print "\n --gff\t\t\t use Kmasker constructed GFF report for annotation";
		print "\n --dbfasta\t\t custom sequences [FASTA] with annotated features in sequence descriptions";
		print "\n --db\t\t pre-calculated blastableDB of nucleotides used for annotation";	
		
		print "\n --hist\t\t create histogram (requires --clist and --occ)";
		print "\n --violin\t\t create violin plot (distribution of k-mer)";
		print "\n --hexplot\t\t create hexagon plot (comparison of two kindex)";
		print "\n --triplot\t\t create tertiary plot (comparison of three kindex)";
		print "\n --occ\t\t provide a Kmasker constructed occ file containing k-mer frequencies";
		print "\n --clist\t\t file containing a list of contig identifier for analysis";	

#		print "\n --stats\t\t\t calculate basic statistics like avegare k-mer frequency per contig etc. (requires --occ)";	
		
		print "\n\n";
		exit();
	}
	

    print "\n Description:\n\t Kmasker is a tool for the automatic detection of repetitive sequence regions.";
    
    print "\n\n Modules:";
	print "\n --build\t\t construction of new index (requires --indexfiles)";
	print "\n --run\t\t\t perform analysis and masking (requires --fasta)";
	print "\n --explore\t\t perform downstream analysis with constructed index and detected repeats";
	
	print "\n\n General options:";
	print "\n --show_repository\t\t shows complete list of global and private k-mer indices";
	print "\n --show_details\t\t\t shows details for a requested kindex";
	print "\n --remove_kindex\t\t remove kindex from repository";
	print "\n --set_private_path\t\t change path to private repository";
	print "\n --set_external_path\t\t change path to external repository [readonly]";
	print "\n --expert_setting_kmasker\t submit individual parameter to Kmasker eg. pctgap,";
	print "\n\t\t\t\t minseed, mingff (see documentation!)";
	print "\n --expert_setting_jelly\t\t submit individual parameter to jellyfish (e.g. on memory usage ";
	print "\n\t\t\t\t for index construction)";
	print "\n --threads\t\t\t set number of threads [4]";
	
	print "\n\n";
	exit();
}

#######
## MAIN

#CHECK settings
if(defined $check_install){
	&check_install();
	exit();
}

#READ global settings
&read_user_config;
&read_repository;


## GET USER INPUT

#THREADS
if(defined $threads_usr){
	$threads = $threads_usr;
}

#index name
if(defined $index_name_usr){
	$index_name = $index_name_usr;
}

#MEM 
if(defined $mem_usr){
	$mem = $mem_usr;
}


########
# STORE GENERAL INFO in HASH info
my %HASH_info 						= ();
$HASH_info{"PATH_kindex_global"}	= $PATH_kindex_global; 
$HASH_info{"PATH_kindex_private"}	= $PATH_kindex_private;
$HASH_info{"PATH_kindex_external"}	= $PATH_kindex_external;
$HASH_info{"path_bin"}				= $path;
$HASH_info{"version KMASKER"}		= $version;
$HASH_info{"temp_path"}				= $temp_path;
$HASH_info{"threads"}		 		= $threads;
$HASH_info{"memory"}			 	= $mem;



#SET private path
if(defined $set_private_path){
	
	print "\n Kmasker will change private path: ".$set_private_path;
	print "\n Please note that data will not be moved to new directories!";
	print "\n This has to be done manually!\n\n";
	
	&set_private_path($set_private_path);
	exit();
}

#SET external path
if(defined $set_external_path){
	
	print "\n Kmasker will change external path: ".$set_external_path;
	print "\n Please note that data will not be moved to new directories!";
	print "\n This has to be done manually!\n\n";
	
	&set_external_path($set_external_path);
	exit();
}

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

#EXPERT SETTINGs Kmasker
if($expert_setting_kmasker ne ""){
	
	my @ARRAY_user_parameter_input = split(";", $expert_setting_kmasker);
	foreach(@ARRAY_user_parameter_input){
		if($_ =~ /:/){
			my @ARRAY_tmp = split(":", $_);
			my $P = $ARRAY_tmp[0];
			if(defined $ARRAY_tmp[1]){
				my $V = $ARRAY_tmp[1];
				
				if($P eq "pctgap"){
					print "\n\n .. changing default parameter for PCTGAP from ".$MK_percent_gapsize." to ". $V." !\n";
					$MK_percent_gapsize = $V;
				}
				
				if($P eq "minseed"){
					print "\n\n .. changing default parameter for MINSEED from ".$MK_min_seed." to ". $V." !\n";
					$MK_min_seed = $V;
				}
				
				if($P eq "mingff"){
					print "\n\n .. changing default parameter for MINGFF from ".$MK_min_gff." to ". $V." !\n";
					$MK_min_gff = $V;
				}
				
			}else{
				print "\n\n Cannot use user input (".$expert_setting_kmasker."). Kmasker is stopped!\n\n\n";
				exit();
			}
		}else{
			print "\n\n Cannot use user input (".$expert_setting_kmasker."). Kmasker is stopped!\n\n\n";
			exit();
		}
	}
}

#genome_size
if(defined $genome_size_usr){
	if($genome_size_usr =~ /^[+-]?\d+$/){
		#is number
		$genome_size = $genome_size_usr;
	}
}

#common_name
if(defined $common_name_usr){
	$common_name = $common_name_usr;
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


#CHECK setting
&check_settings;

if(defined $build){
	#USE BUILD MODULE
	
	if(defined $set_global){
		if(exists $HASH_repository_kindex{$set_global}){
			&set_kindex_global($set_global, \%HASH_repository_kindex, $PATH_kindex_private, $PATH_kindex_global);			
		}else{
			print "\n WARNING: requested kindex '".$set_global."' does not exist. Kmasker was stopped!\n\n";
		}
		exit();			
	}else{
		
		########
		#STORE BUILD INFO in HASH info
		my $input 							= join(" ", sort { $a cmp $b } @seq_usr);
		$HASH_info{"user name"}				= $user_name;
		$HASH_info{"seq"} 					= $input;		
		#REQUIRED
		$HASH_info{"k-mer"}					= $k 			if(defined $k);
		$HASH_info{"genome size"}			= $genome_size	if(defined $genome_size);
		$HASH_info{"kindex name"}			= $index_name	if(defined $index_name);
		$HASH_info{"common name"}			= $common_name  if(defined $common_name);
		#ADDITIONAL
		$HASH_info{"expert setting kmasker"}= $expert_setting_kmasker; 
		$HASH_info{"expert setting jelly"}	= $expert_setting_jelly; 
		$HASH_info{"version KMASKER"}		= $version;
		$HASH_info{"version BUILD"} 		= "";
		$HASH_info{"status"}				= "";
		$HASH_info{"scientific name"}		= "";
		$HASH_info{"sequence type"}			= "";
		$HASH_info{"general notes"}			= "";
		$HASH_info{"type"}					= "";
		$HASH_info{"sequencing depth"}		= "";
		
		#CONSTRUCT
		if(defined $input){
			&build_kindex_jelly(\%HASH_info, $build_config, \%HASH_repository_kindex, \%HASH_path); 
		}
		
		#CLEAN
		&clean_repository_directory(\%HASH_info, \%HASH_repository_kindex);			
	}	
		
	#QUIT
	print "\n - Thanks for using Kmasker! -\n\n";
	exit();
}

if(defined $run){
	#USE RUN MODULE
	
	########
	#STORE RUN INFO in HASH info
	$HASH_info{"user_name"}				= $user_name;
	$HASH_info{"kindex name"}			= $kindex;
	$HASH_info{"rept"}					= $repeat_threshold;
	$HASH_info{"min_length"}			= $length_threshold; 
	$HASH_info{"MK_percent_gapsize"}	= $MK_percent_gapsize;
	$HASH_info{"MK_min_seed"}			= $MK_min_seed;
	$HASH_info{"MK_min_gff"}			= $MK_min_gff;
	$HASH_info{"expert setting kmasker"}= $expert_setting_kmasker; 
	$HASH_info{"expert setting jelly"}	= $expert_setting_jelly; 
	$HASH_info{"version KMASKER"}		= $version;
	$HASH_info{"version BUILD"} 		= "";

	
	if(defined $kindex){
	#single kindex		
	
		#READ repository.info
		my $FILE_repository_info = "";
		if(exists $HASH_repository_kindex{$kindex}){
			my @ARRAY_repository	= split("\t", $HASH_repository_kindex{$kindex});
			$FILE_repository_info 	= $ARRAY_repository[4]."repository_".$kindex.".info";
		}else{
			print "\n .. Kmasker was stopped. Info for kindex does not exists!\n";
			exit ();
		}
	
		my $href_info 	= &read_config($FILE_repository_info, \%HASH_info, \%HASH_repository_kindex, "run");
		%HASH_info		= %{$href_info};
		
		#START RUN			
		&run_kmasker_SK($fasta, $kindex, \%HASH_info, \%HASH_repository_kindex);

	}elsif(scalar(@multi_kindex > 1)){
	#multiple kindex
	
		my @ARRAY_HASH_info_aref = ();
		
		for(my $ki=0;$ki<scalar(@multi_kindex);$ki++){
			#READ repository.info
			my $kindex_K = $multi_kindex[$ki];
			my $FILE_repository_info = "";
			if(exists $HASH_repository_kindex{$kindex_K}){
				my %HASH_info_Kx 				= %HASH_info;
				$HASH_info_Kx{"kindex name"}	= $kindex_K;
				my @ARRAY_repository			= split("\t", $HASH_repository_kindex{$kindex_K});
				$HASH_info_Kx{"k-mer"}			= $ARRAY_repository[1];
				$FILE_repository_info 			= $ARRAY_repository[4]."repository_".$kindex_K.".info";				
				
				#get all infos for repository
				my $href_info_Kx 	= &read_config($FILE_repository_info, \%HASH_info_Kx, \%HASH_repository_kindex, "run");
				$ARRAY_HASH_info_aref[$ki] 		= $href_info_Kx;
			}else{
				print "\n .. Kmasker was stopped. Info for kindex (".$kindex_K.") does not exists!\n";
				exit ();
			}			
		}
	
		&run_kmasker_MK($fasta, \@multi_kindex, \@ARRAY_HASH_info_aref, \%HASH_repository_kindex);
	}
	
	#QUIT
	print "\n - Thanks for using Kmasker! -\n\n";
	exit();
}

if(defined $explore){
	#USE EXPLORE MODULE
	
	########
	#STORE EXPLORE INFO in HASH info
	$HASH_info{"user name"}				= $user_name;
	$HASH_info{"kindex"}				= $kindex			if(defined $kindex);
	$HASH_info{"multi kindex"}			= \@multi_kindex 	if(defined $multi_kindex[0]);
	#ADDITIONAL
	$HASH_info{"expert setting kmasker"}= $expert_setting_kmasker; 
	$HASH_info{"expert setting jelly"}	= $expert_setting_jelly; 
	$HASH_info{"expert setting blast"}	= $expert_setting_blast if(defined $expert_setting_blast); 
	$HASH_info{"version KMASKER"}		= $version;
	$HASH_info{"version BUILD"} 		= "";	
	
	
	#ANNOTATION
	if(defined $custom_annotate){
		# EXPLORE ANNOTATE by GFF
		
		#START annotation of sequence features, provided in GFF with custome FASTA sequence or blastableDB
		my $check_settings = 1;
		$check_settings = 0 if(!defined $gff);
		$check_settings = 0 if(!defined $fasta);
				
		if($check_settings == 1){
			my %HASH_db;
    		$HASH_db{"dbfasta"} = $dbfasta if(defined $dbfasta);
    		$HASH_db{"db"} 		= $dbfasta if(defined $blastableDB);
			&custom_annotation($fasta, $gff, \%HASH_db, \%HASH_info);
		}else{
			print "\n WARNING: Required parameter are missing. Kmasker was stopped.\n\n";
			exit;
		}		
	}
	
	my $visualisation;
	$visualisation = 1 if(defined $hist);
	$visualisation = 1 if(defined $violin);
	$visualisation = 1 if(defined $hexplot);
	$visualisation = 1 if(defined $triplot);
		
	#HISTOGRAM
	if(defined $visualisation){
		# EXPLORE VIZUALISATIONS
	
		if(defined $occ){
		# explore requires an OCC file
			
			my $missing_parameter = "";
		
			if(defined $hist){
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
			print "\n ERROR: no occ provided. For this 'Kmasker --explore' subfunction an occ file is required!\n\n";
		}
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
	&remove_kindex($remove_kindex,\%HASH_repository_kindex);
	exit();
}				

##END MAIN


## subroutine
#
sub read_user_config(){
	$user_name 			= `whoami`;
	$user_name			=~ s/\n//g;
	my $gconf 			= $path."/kmasker.config";
	my $uconf 			= $ENV{"HOME"}."/.kmasker_user.config";
	
	if(-e $gconf){
		#LOAD global info
		my $gCFG = new IO::File($gconf, "r") or die "\n unable to read user config $!";	
		
		while(<$gCFG>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_tmp = split("=", $line);
			$PATH_kindex_global	= $ARRAY_tmp[1] if($ARRAY_tmp[0] eq "PATH_kindex_global");
			$PATH_kindex_global .= "/" if($PATH_kindex_global !~ /\/$/);
			
			#GLOBAL
			#PRIVATE
			if(-d $PATH_kindex_global){
				#directory exists - do nothing
			}else{
				#directory has to be created
				system("mkdir "."\"".$PATH_kindex_global."\"");
			}
			
			#READ external tool path
			#JELLYFISH
			if($line =~ /^jellyfish=/){
				my @ARRAY_tmp = split("=", $line);
				if(!defined $ARRAY_tmp[1]){
					system("which jellyfish >/dev/null 2>&1 || { echo >&2 \"Kmasker requires jellyfish but it's not installed! Kmasker process stopped.\"; exit 1; \}");
					$HASH_path{"jellyfish"} = `which jellyfish`;
					$HASH_path{"jellyfish"} =~ s/\n//;
				}else{
					$HASH_path{"jellyfish"} = $ARRAY_tmp[1];
					system("which ".$HASH_path{"jellyfish"}." >/dev/null 2>&1 || { echo >&2 \"Kmasker requires jellyfish but it's not installed!  Kmasker process stopped.\"; exit 1; \}");
				}
				print "\n jellyfish=".$HASH_path{"jellyfish"}."\n" if(defined $verbose);
			}
			
			#FASTQ-STATs
			if($line =~ /^fastq-stats=/){
				my @ARRAY_tmp = split("=", $line);
				if(!defined $ARRAY_tmp[1]){
					system("which fastq-stats >/dev/null 2>&1 || { echo >&2 \"Kmasker requires fastq-stats but it's not installed! Kmasker process stopped.\"; exit 1; \}");
					$HASH_path{"fastq-stats"} = `which fastq-stats`;
					$HASH_path{"fastq-stats"} =~ s/\n//;
				}else{
					$HASH_path{"fastq-stats"} = $ARRAY_tmp[1];
					system("which ".$HASH_path{"fastq-stats"}." >/dev/null 2>&1 || { echo >&2 \"Kmasker requires fastq-stats but it's not installed! Kmasker process stopped.\"; exit 1; \}");
				}
				print "\n fastq-stats=".$HASH_path{"fastq-stats"}."\n" if(defined $verbose);
			}
			
			#GFFREAD
			if($line =~ /^gffread=/){
				my @ARRAY_tmp = split("=", $line);
				if(!defined $ARRAY_tmp[1]){
					system("which gffread >/dev/null 2>&1 || { echo >&2 \"Kmasker requires gffread but it's not installed! Kmasker process stopped.\"; exit 1; \}");
					$HASH_path{"gffread"} = `which gffread`;
					$HASH_path{"gffread"} =~ s/\n//;
				}else{
					$HASH_path{"gffread"} = $ARRAY_tmp[1];
					system("which ".$HASH_path{"gffread"}." >/dev/null 2>&1 || { echo >&2 \"Kmasker requires gffread but it's not installed! Kmasker process stopped.\"; exit 1; \}");
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
			my @ARRAY_tmp = split("=", $line);
			$PATH_kindex_private= $ARRAY_tmp[1] if($ARRAY_tmp[0] eq "PATH_kindex_private");
			$PATH_kindex_private.= "/" if($PATH_kindex_private !~ /\/$/);
			
			#PRIVATE
			if(-d $PATH_kindex_private){
				#directory exists - do nothing
			}else{
				#directory has to be created
				system("mkdir "."\"".$PATH_kindex_private."\"");
			}			
			
			#EXTERNAL
			if(($ARRAY_tmp[0] eq "PATH_kindex_external")&&(defined $ARRAY_tmp[1])){
				$PATH_kindex_external = $ARRAY_tmp[1];
				$PATH_kindex_external.= "/" if($PATH_kindex_external !~ /\/$/);
				if($PATH_kindex_external !~ /^\/$/){
					if(-d $PATH_kindex_external){
						#directory exists - do nothing
					}else{
						#directory has to be created
						system("mkdir "."\"".$PATH_kindex_external."\"");
					}	
				}		
			}					
		}
	}else{
		#SETUP user
		&initiate_user();
	}
}


## subroutine
#
sub read_repository(){
	
	#PRIVATE
	opendir( my $DIR_P, $PATH_kindex_private ) or die "Can not open \'$PATH_kindex_private\' (private path)\n";
	my $status 				= "";
	my $common_name_this 	= "";
	my $file_name			= "";
	my $kmer				= "";
	while ( $file_name = readdir $DIR_P ) {
		$status 			= "private";
		$common_name_this 	= "";		
		if($file_name =~ /^repository_/){
			$file_name =~ s/repository_//;
			my @ARRAY_name 	= split(/\./, $file_name);
			my $kindex_id 			= $ARRAY_name[0]; 
			my $BUILD_file 	= new IO::File($PATH_kindex_private."repository_".$kindex_id.".info", "r") or print " ... could not read repository info for $kindex_id : $!\n";
			if(-e $PATH_kindex_private."repository_".$kindex_id.".info"){;
				while(<$BUILD_file>){
					if($_ =~ /^common name/){
						$common_name_this = +(split("\t", $_))[1];
						$common_name_this =~ s/\n//;
					}					
					if($_ =~ /^k-mer/){
						$kmer = +(split("\t", $_))[1];
						$kmer =~ s/\n//;
					}					
				}
				#integrate into HASH
				$HASH_repository_kindex{$kindex_id} = $kindex_id."\t".$kmer."\t".$common_name_this."\t".$status."\t".$PATH_kindex_private;
			}
		}
	}
	close $DIR_P;
	
	#EXTERNAL
	if($PATH_kindex_external !~ //){
	# EXTERNAL PATH is set
		opendir( my $DIR_E, $PATH_kindex_external ) or die "Can not open \'$PATH_kindex_external\' (external path)\n";
		$common_name_this 	= "";
		$file_name			= "";
		$kmer				= "";
		while ( $file_name = readdir $DIR_E ) {
			$status 			= "external";
			$common_name_this 	= "";		
			if($file_name =~ /^repository_/){
				$file_name =~ s/repository_//;
				my @ARRAY_name 	= split(/\./, $file_name);
				my $kindex_id 			= $ARRAY_name[0]; 
				my $BUILD_file 	= new IO::File($PATH_kindex_external."repository_".$kindex_id.".info", "r") or print " ... could not read repository info for $kindex_id : $!\n";
				if(-e $PATH_kindex_external."repository_".$kindex_id.".info"){;
					while(<$BUILD_file>){
						if($_ =~ /^common name/){
							$common_name_this = +(split("\t", $_))[1];
							$common_name_this =~ s/\n//;
						}					
						if($_ =~ /^k-mer/){
							$kmer = +(split("\t", $_))[1];
							$kmer =~ s/\n//;
						}					
					}
					#integrate into HASH
					$HASH_repository_kindex{$kindex_id} = $kindex_id."\t".$kmer."\t".$common_name_this."\t".$status."\t".$PATH_kindex_external;
				}
			}
		}
		close $DIR_E;
	}
	
	#GLOBAL
	opendir( my $DIR_G, $PATH_kindex_global ) or die "Can not open $PATH_kindex_private\n"; ;
	$common_name_this = "";
	while ( $file_name = readdir $DIR_G ) {
		$status 			= "global";
		$common_name_this 	= "";	
		if($file_name =~ /^repository_/){
			$file_name =~ s/repository_//;
			my @ARRAY_name 	= split(/\./, $file_name);
			my $kindex_id 			= $ARRAY_name[0]; 
			my $BUILD_file 	= new IO::File($PATH_kindex_global."repository_".$kindex_id.".info", "r") or print " ... could not read repository info for $kindex_id : $!\n";
			if(-e $PATH_kindex_global."repository_".$kindex_id.".info"){;
				while(<$BUILD_file>){
					if($_ =~ /^common name/){
						$common_name_this = +(split("\t", $_))[1];
						$common_name_this =~ s/\n//;
					}	
					if($_ =~ /^k-mer/){
						$kmer = +(split("\t", $_))[1];
						$kmer =~ s/\n//;
					}	
				}
				#integrate into HASH
				$HASH_repository_kindex{$kindex_id} = $kindex_id."\t".$kmer."\t".$common_name_this."\t".$status."\t".$PATH_kindex_global;
			}		
		}
	}
	
	close $DIR_G;
}


## subroutine
#
sub show_repository(){	
	print "\n\nREPOSITORY of available kindex structures:\n";		
	#PRINT
	foreach my $kindex_this (keys %HASH_repository_kindex){
		my @ARRAY_repository	= split("\t", $HASH_repository_kindex{$kindex_this});
		#formatted print
		print "\n\t";
		printf "%-14s", $kindex_this;				
		print "\t".$ARRAY_repository[1]."\t";
		printf "%-14s", $ARRAY_repository[2];
		print "\t".$ARRAY_repository[3];
	}
	print "\n\n";
}

## subroutine
#
sub show_details_for_kindex(){
	my $kindex = $_[0];
	if(exists $HASH_repository_kindex{$kindex}){
		my @ARRAY_repository	= split("\t", $HASH_repository_kindex{$kindex});
		my $BUILD_file 		= new IO::File($ARRAY_repository[4]."repository_".$kindex.".info") or die " ... can not read/find repository_".$kindex.".info file in ".$ARRAY_repository[4]." details : $!\n\n";
		print "\n  KINDEX details for ".$kindex.": \n";
		while(<$BUILD_file>){
			print "\t".$_;
		}
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
	my $uconf 			= $ENV{"HOME"}."/.kmasker_user.config";

	if(-e $uconf){
		#USER already exists, do nothing
	}else{
		#SETUP user conf
		my $USER_CONF 	= new IO::File($uconf, "w") or die "could not write user repository : $!\n";
		print $USER_CONF "PATH_kindex_private=".$ENV{"HOME"}."/KINDEX/";
		print $USER_CONF "\nPATH_kindex_external=";
		close $USER_CONF;	
		
		#SHOW INFO
		print "\n PLEASE NOTE: \n You are writing all large data structures to your home directory [default].";
		print "\n It is recommended to modify the path for 'PATH_kindex_private'.\n";
		print "\n Use the following command: 'Kmasker --set_private_path enter/your/path'\n\n";	
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
	if(defined $explore){
		$module_count++; 
	}
	
	#MULTIPLE 
	if($module_count > 1){
		print "\n Kmasker was stopped. Multiple modules (build, run or explore were used!\n";
		exit(0);
	}	
}


## subroutine
#
sub check_install(){

	$user_name 			= `whoami`;
	$user_name			=~ s/\n//g;
	my $gconf 			= $path."/kmasker.config";
	
	#PERMISSION - calling this procedure is only be possible for directory owner (who installed Kmasker)
	my $fp			 =  $path."/kmasker.config";
	my $installed_by = `stat -c "%U" $fp`;
	$installed_by =~ s/\n//;
	if($installed_by ne $user_name){
		print "\n Your user rights are not sufficient to call that procedure. Call is permitted.\n";
		print "I=(".$installed_by.") U=(".$user_name.")\n";
		exit();	
	}
	
	#REQUIREMENTs
	my %HASH_requirments	= ("PATH_kindex_global" => $ENV{"HOME"}."/KINDEX/",
								"jellyfish" => "",
								"fastq-stats" => "",
								"gffread" => "");
			
	#SET default path if tool is detected
	foreach my $tool (keys %HASH_requirments){
		if($HASH_requirments{$tool} eq ""){
			$HASH_requirments{$tool} = `which $tool`;
			$HASH_requirments{$tool} =~ s/\n//;
			print "\n DEFAULT (".$tool.")= ".$HASH_requirments{$tool};
		}
	}
	
	#GLOBAL
	if(-e $gconf){
		#LOAD global info
		my $gCFG_old 	= new IO::File($gconf, "r") or die "\n unable to read user config $!";	
		my $gCFG 		= new IO::File($gconf.".tmp", "w") or die "\n unable to update user config $!";
		
		my %HASH_provided = ();
		while(<$gCFG_old>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_tmp = split("=", $line);
			
			#PATH
			if($line =~ /^PATH_kindex_global/){
				if(defined $ARRAY_tmp[1]){
					$HASH_requirments{"PATH_kindex_global"} = $ARRAY_tmp[1];
					$HASH_requirments{"PATH_kindex_global"}	.= "/" if($HASH_requirments{"PATH_kindex_global"} !~ /\/$/ );
				}
			}
			
			$HASH_provided{"jellyfish"} 	= $line if($line =~ /^jellyfish=/);
			$HASH_provided{"fastq-stats"} 	= $line if($line =~ /^fastq-stats=/);
			$HASH_provided{"gffread"}		= $line if($line =~ /^gffread=/);			
		}
		
		
		#CHECK tool requirments
		#JELLYFISH
		$HASH_requirments{"jellyfish"} = &check_routine_for_requirement("jellyfish", $HASH_provided{"jellyfish"}, $HASH_requirments{"jellyfish"});
		
		#FASTQ-STATs
		$HASH_requirments{"fastq-stats"} = &check_routine_for_requirement("fastq-stats", $HASH_provided{"fastq-stats"}, $HASH_requirments{"fastq-stats"});
		
		#GFFREAD
		$HASH_requirments{"gffread"} = &check_routine_for_requirement("gffread", $HASH_provided{"gffread"}, $HASH_requirments{"gffread"});
				
		#WRITE
		print $gCFG "#path requirements\n";
		print $gCFG "PATH_kindex_global=".$HASH_requirments{"PATH_kindex_global"}."\n\n";
		print $gCFG "#external tool requirements\n";
		foreach my $required (keys %HASH_requirments){
			if($required !~ /^PATH_kindex/){
				system("which $HASH_requirments{$required} >/dev/null 2>&1 || { echo >&2 \"Kmasker requires $required but it's not installed or path is missing! Kmasker process stopped.\"; exit 1; \}");
				print $gCFG $required."=".$HASH_requirments{$required}."\n";
				print "\n write ".$required." --> ".$HASH_requirments{$required};
			}			
		}
		
		if($HASH_requirments{"PATH_kindex_global"} eq $ENV{"HOME"}."/KINDEX/"){					
			print "\n PATH_kindex_global=".$ENV{"HOME"}."/KINDEX/\n\n";
			print "\n The global path variable is not defined and was automatically set to your home directory."; 
			print "\n Its recommended to use another directoty because large data volumes will be produced and stored in that directory.";
			print "\n Please edit the global path variable for storage of constructed KINDEX by using the following command";
			print "\n\n  Kmasker --set_global_path /global_path/\n\n"; 
		}
		
		$gCFG_old->close();
		$gCFG->close();
		system("mv ".$gconf.".tmp ".$gconf)	
	}
}

sub check_routine_for_requirement(){
	my $requirement = $_[0];
	my $line 		= $_[1];
	my $default		= $_[2];
	
	if($line =~ /^$requirement=/){
		my @ARRAY_tmp = split("=", $line);
		if(defined $ARRAY_tmp[1]){
			#CHECK if path is correct
			my $path_given 	= $ARRAY_tmp[1];
			my $path_check 	= `which $path_given`;
			$path_check		=~ s/\n$//;
			if($path_check eq ""){
				print "\n ... provided path for ".$requirement." seems to be wrong! Trying to detect path automatically\n";			
			}else{
				$default = $path_check;
			}					
		}
	}	
	return $default;
}


1;
