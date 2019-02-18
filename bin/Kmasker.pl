#!/usr/bin/perl -w
use strict;
use warnings;
use IO::File;
use Getopt::Long;
#use Digest::MD5 qw(md5 md5_hex md5_base64);
use File::stat;
use File::Copy qw(move);
#setup package directory
use File::Basename;
use POSIX qw(strftime);
use Cwd  qw(abs_path getcwd);
use lib dirname(dirname abs_path $0) . '/lib';

#include packages
use kmasker::kmasker_build qw(build_kindex_jelly remove_kindex set_kindex_external set_private_path set_external_path show_path_infos clean_repository_directory read_config);
use kmasker::kmasker_run qw(run_kmasker_SK run_kmasker_MK run_krispr show_version_PM_run);
use kmasker::kmasker_explore qw(plot_histogram_raw plot_histogram_mean custom_annotation report_statistics plot_maker plot_maker_direct plot_barplot);
use kmasker::functions;

my $version 	= "0.0.35 rc190218";
my $path 		= dirname abs_path $0;		
my $indexfile;
my $PID = $kmasker::functions::PID;
my $uPID; 
my $LONG_PID;

#MODULES
#BUILD
my $build;
my $run;
my $explore;
my $repositories;
my $build_config;
my $make_config;
my $seq_type 			= "reads";
my $seq_type_usr;
my $genome_size;
my $genome_size_usr;
my $scientific_name 	= "";
my $scientific_name_usr;
my $common_name 		= "";
my $common_name_usr ;
my $index_name;
my $index_name_usr;
my $threads_usr;
my $threads = 4;
my $size_usr;
my $size			= 24;
my $foldchange		= 10;
my $foldchange_usr;
my $PATH_kindex_private = "";
my $PATH_kindex_external= "";

#RUN
my $fasta;
my $krispr;
my $fish;
my $compare;
my $bed;
my $k 					= 21;
my $k_usr;
my $mismatch			= 0;
my $mismatch_usr;
my $tool_jellyfish;
my @seq_usr;
my $length_threshold	= 100;
my $length_threshold_usr;
my $repeat_threshold	= 10;
my $repeat_threshod_usr;
my $tolerant_length_threshold_usr;
my $tolerant_length_threshold = 0;

# set the following parameter using '--expert_setting_kmasker' or provide them in file using '--config_kmasker'
my $MK_percent_gapsize	= 10;	#Default	MK_percent_gapsize (N) is the paramater allowing N percent of a region to be not repetitive when merging regions. Used in SK and MK.
my $MK_min_seed			= 5;	#Default	MK_min_seed (N) is the minimal length (bp) of a featuere to start in the process of merging regions for GFF
my $MK_min_gff			= 25;   #Default	MK_min_gff (N) is the minimal length (bp) of a featue to be reported in GFF

#EXPLORE
my $gff;
my $list;
my $file;
my $cfile;
my $occ;
my @OCClist;
my $stats;
my $cstats;	
my $custom_annotate;
my $blastableDB;
my $dbfasta;
my $dynamic;
my $force;
my $xtract;
#visualisation
my $hist;
my $histm;
my $violin;
my $boxplot;
my $hexplot;
my $barplot;
my $sws;
my $max;
my $bin;
#my $log;
#GENERAL parameter
my $help;
my $out;
my $keep_temporary_files;
my $show_kindex_repository;
my $show_details_for_kindex;
my $show_path;
my $show_version;
my $set_private_path;
my $set_external_path;
my $check_install;
my $remove_kindex;
my $set_global;
my $user_name = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
my $verbose;
my $temp_path			= "./temp/";
my $feature 			= "KRC";

#CONFIGURATION parameter
my $expert_setting_kmasker 	= "",
my $expert_setting_jelly 	= ""; 
my $expert_setting_blast	= "";
my $expert_config_kmasker;
my $expert_config_jelly;
my $expert_config_blast;
my $gconf = dirname($path)."/config/kmasker.config";
my $uconf = $ENV{"HOME"}."/.kmasker_user.config";
my @ARGV_B = @ARGV;


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
							"sn=s"				=> \$scientific_name_usr,
							"cn=s"				=> \$common_name_usr,
							"in=s"				=> \$index_name_usr,
							"config=s"			=> \$build_config,
							"as"				=> \$seq_type_usr,
							"make_config"		=> \$make_config,							
							
							#RUN
							"fasta=s"			=> \$fasta,	
							"krispr"			=> \$krispr,
							"fish"				=> \$fish,
							"compare"			=> \$compare,
							"kindex=s{1,}"		=> \@multi_kindex,
							"rept=i"			=> \$repeat_threshod_usr,
							"minl=i"			=> \$length_threshold_usr,
							"m=i"				=> \$mismatch,							
											
							#EXPLORE
							"annotate"			=> \$custom_annotate,
							"gff=s"				=> \$gff,
							"feature=s"			=> \$feature,
							"db=s"				=> \$blastableDB,
							"dbfasta=s"			=> \$dbfasta,
							"barplot"			=> \$barplot,
							"boxplot"			=> \$boxplot,
							"hexplot"			=> \$hexplot,		
							"hist"				=> \$hist,
							"histm"				=> \$histm,
							"violin"			=> \$violin,
							"cfile=s"			=> \$cfile,
							"file=s"			=> \$file,
							"list=s"			=> \$list,
							"occ=s{1,}"			=> \@OCClist,
							"stats"				=> \$stats,	
							"cstats"			=> \$cstats,	
							"dynamic"			=> \$dynamic,
							"window"			=> \$sws,	
							"max=i"				=> \$max,	
							"bin=i"				=> \$bin,
#							"log"				=> \$log,	
							"xtract"			=> \$xtract,				
							
							#GLOBAL
							"show_repository"			=> \$show_kindex_repository,
							"show_details=s"			=> \$show_details_for_kindex,
							"show_path"					=> \$show_path,
							"remove_kindex=s"			=> \$remove_kindex,
							"set_private_path=s"		=> \$set_private_path,
							"set_external_path=s"		=> \$set_external_path,							
							"check_install"				=> \$check_install,	
							"threads=i"					=> \$threads_usr,
							"out=s"						=> \$out,
							
							#configuration
							"expert_setting_jelly=s"	=> \$expert_setting_jelly,
							"expert_setting_kmasker=s"	=> \$expert_setting_kmasker,
							"expert_setting_blast=s"	=> \$expert_setting_blast,
							"config_jelly=s"			=> \$expert_config_jelly,
							"config_kmasker=s"			=> \$expert_config_kmasker,
							"config_blast=s"			=> \$expert_config_blast,
							"force"						=> \$force,
							"bed"						=> \$bed,
							"global_conf"				=> \$gconf,
							"user_conf"					=> \$uconf,
							#Houskeeping
							"keep_tmp"					=> \$keep_temporary_files,
							"verbose"					=> \$verbose,
							"show_version"				=> \$show_version,
							"long_pid"					=> \$LONG_PID,
							"pid=s"						=> \$uPID,
							"help"						=> \$help										
						);
						
#READ global settings
&read_user_config;
&read_repository;
#&create_PID;



#######
## MAIN


## GET USER INPUT

#THREADS
if(defined $threads_usr){
	$threads = $threads_usr;
}

#index name
if(defined $index_name_usr){
	$index_name = $index_name_usr;
}

#sequence type 
if(defined $seq_type_usr){
	$seq_type = "assembly";
}
if(defined $LONG_PID){
	kmasker::functions::set_long_PID();
	$PID = $kmasker::functions::PID;
	$kmasker::kmasker_run::log="log_run_" . $kmasker::functions::PID . ".txt";
	$kmasker::kmasker_explore::log="log_explore_" . $kmasker::functions::PID . ".txt";
	$kmasker::kmasker_explore::PID=$kmasker::functions::PID;
}
if(defined $uPID) {
	$PID = $uPID;
	$kmasker::functions::PID = $uPID;
	$kmasker::kmasker_run::log="log_run_" . $uPID . ".txt";
	$kmasker::kmasker_explore::log="log_explore_" . $uPID . ".txt";
	$kmasker::kmasker_explore::PID=$uPID;
}

#############
# CALLING general commands and tasks
&intro_call();					
			
			
#############
# CALLING HELP if requested						
$help = 1 if((!defined $build)&&(!defined $run)&&(!defined $explore));

if(defined $help){
	&help();	
}

########
# STORE GENERAL INFO in HASH info
my %HASH_info 						= ();
$HASH_info{"PATH_kindex_private"}	= $PATH_kindex_private;
$HASH_info{"PATH_kindex_external"}	= $PATH_kindex_external;
$HASH_info{"path_bin"}				= $path;
$HASH_info{"version KMASKER"}		= $version;
$HASH_info{"temp_path"}				= $temp_path;
$HASH_info{"threads"}		 		= $threads;
$HASH_info{"verbose"}				= $verbose if(defined $verbose);
$HASH_info{"MK_percent_gapsize"}	= $MK_percent_gapsize;
$HASH_info{"MK_min_seed"}			= $MK_min_seed;
$HASH_info{"MK_min_gff"}			= $MK_min_gff;
$HASH_info{"bed"}					= 0;
$HASH_info{"bed"}					= 1 if(defined $bed);
$HASH_info{"fold-change"}			= $foldchange;
$HASH_info{"expert setting kmasker"}= $expert_setting_kmasker; 
$HASH_info{"expert setting jelly"}	= $expert_setting_jelly;
$HASH_info{"expert setting blast"}	= $expert_setting_blast;
$HASH_info{"size"}					= $size;
$HASH_info{"krisp mismatch"}		= $mismatch;
$HASH_info{"PID"}					= $PID;
$HASH_info{"LONG_PID"}				= $LONG_PID;


########
# STORE INFO about DB in HASH_db
my %HASH_db 						= ();


#USER specification
#kindex
if(scalar(@multi_kindex) == 1){
	if(exists $HASH_repository_kindex{$multi_kindex[0]}){
		$kindex = $multi_kindex[0];
	}else{
		print "\n ERROR: defined kindex ('".$multi_kindex[0]."') does not exist!\n\n";
		exit();
	}	
}

#EXPERT SETTING from configuration file
&read_configuration_file("blast", $expert_config_blast) if(defined $expert_config_blast);	
&read_configuration_file("jelly", $expert_config_jelly) if(defined $expert_config_jelly);	
&read_configuration_file("kmasker", $expert_config_kmasker) if(defined $expert_config_kmasker);	

#EXPERT SETTING from command line Kmasker
&use_expert_settings("kmasker", $expert_setting_kmasker) if($expert_setting_kmasker ne "");	
&use_expert_settings("jelly", $expert_setting_jelly) if($expert_setting_jelly ne "");	
&use_expert_settings("blast", $expert_setting_blast) if($expert_setting_blast ne "");	

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

#scientific name
if(defined $scientific_name_usr){
	$scientific_name = $scientific_name_usr;
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


#######################
###
### BUILD SECTION
###
if(defined $build){
	#USE BUILD MODULE
	
	if(defined $index_name){
		if(exists $HASH_repository_kindex{$index_name}){
			print "\n WARNING: KINDEX ".$index_name." already exists! Kmasker was stopped\n\n";
			exit();
		}	
			
		if(length($index_name) >=25){
			print "\n WARNING: KINDEX name ".$index_name." is to long. Please use index name with less than 25 letters!\n\n";
			exit();
		}	
	}		
		
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
	$HASH_info{"version KMASKER"}		= $version;
	$HASH_info{"version BUILD"} 		= "";
	$HASH_info{"status"}				= "";
	$HASH_info{"scientific name"}		= $scientific_name;
	$HASH_info{"sequence type"}			= $seq_type;
	$HASH_info{"general notes"}			= "";
	$HASH_info{"type"}					= "";
	$HASH_info{"sequencing depth"}		= "";
		
	print "\n starting build module ... \n";	
		
	#CONSTRUCT
	if(defined $input){
		&build_kindex_jelly(\%HASH_info, $build_config, \%HASH_repository_kindex, \%HASH_path); 
	}
		
	#CLEAN
	&clean_repository_directory(\%HASH_info, \%HASH_repository_kindex);			

		
	#QUIT
	print "\n - Thanks for using Kmasker! -\n\n";
	exit();
}


#######################
###
### RUN SECTION
###
if(defined $run){
	#USE RUN MODULE
	
	########
	#STORE RUN INFO in HASH info
	$HASH_info{"user_name"}				= $user_name;
	$HASH_info{"kindex name"}			= $kindex;
	$HASH_info{"rept"}					= $repeat_threshold;
	$HASH_info{"minl"}					= $length_threshold; 
	$HASH_info{"version KMASKER"}		= $version;
	$HASH_info{"version BUILD"} 		= "";
	
	print "\n starting run module ... session id: $PID\n";
	
	#CRISPR design modul
	if(defined $krispr){
		if(exists $HASH_repository_kindex{$kindex}){
			&run_krispr($kindex, $fasta, \%HASH_repository_kindex, \%HASH_info);
		}
		exit();	
	}
	
	
	#FISH design modul
	if(defined $fish){
		
		if(scalar(@multi_kindex) == 1){
			#single kindex	
			
			if(exists $HASH_repository_kindex{$kindex}){
				
				if(!defined $repeat_threshod_usr){
					$repeat_threshold = 10;
				}				
				if(defined $length_threshold_usr){
					if($length_threshold < 2500){
						print "\n WARNING: Provided length is too small for proper FISH candidate sequence. Setting '--minl' to 2500 bp!\n\n";
						$length_threshold = 2500;
					}
				}else{
						$length_threshold = 2500;
				}
				
				my $this_setting = "--rept=".$repeat_threshold.";--minl=".$length_threshold;
				$HASH_info{"rept"} = $repeat_threshold;
				$HASH_info{"minl"} = $length_threshold;
				&use_expert_settings("kmasker", $this_setting);
				
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
				
				#CHECK FISH results
				my $call = "grep \">\" -c KMASKER_filtered_regions_KDX_${kindex}_${PID}.fasta";
				my $this = `$call`;

				
				#CALL PLOT routines
				
				my $occ_kmer_counts = "KMASKER_kmer_counts_KDX_${kindex}_${PID}.occ";
				if($verbose) {
					system("Kmasker --explore --verbose --hist --occ ".$occ_kmer_counts);
					system("Kmasker --explore --verbose --histm --occ ".$occ_kmer_counts);	
				}
				else{
					system("Kmasker --explore --hist --occ ".$occ_kmer_counts);
					system("Kmasker --explore --histm --occ ".$occ_kmer_counts);	
				}	
				
				$this =~ s/\n//;
				if($this > 0){
					print "\n .. Kmasker detected ".$this." candidates for FISH !!\n\n";
					move("mv KMASKER_filtered_regions_KDX_${kindex}_${PID}.fasta", "KMASKER_filtered_fish_KDX_".$kindex."_".$PID.".fasta");	
					
				}else{
					print "\n .. no candidates for FISH detected!!\n\n";
					#system("rm KMASKER_filtered_regions_*".$fasta);
				}
				
				
											
			}
		}else{
			print "\n WARNING: Calling '--fish' settings is only permitted with single KINDEX.";
			print "\n          Kmasker was stopped!\n\n"; 
		}
		exit();	
	}

	
	# K-mer analysis
	if(scalar(@multi_kindex) == 1){
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

	}else{
	#multiple kindex
	
		my @ARRAY_HASH_info_aref = ();
		
		#QUICK check existence if kindex
		for(my $ki=0;$ki<scalar(@multi_kindex);$ki++){
			#READ repository.info
			my $kindex_K = $multi_kindex[$ki];
			if(!exists $HASH_repository_kindex{$kindex_K}){
				print "\n .. Kmasker was stopped. Info for kindex (".$kindex_K.") does not exists!\n";
				exit(); 
			}
		}		
		
		for(my $ki=0;$ki<scalar(@multi_kindex);$ki++){
			#READ repository.info
			my $kindex_K = $multi_kindex[$ki];
			my $FILE_repository_info = "";
			my %HASH_info_Kx 				= %HASH_info;
			$HASH_info_Kx{"kindex name"}	= $kindex_K;
			my @ARRAY_repository			= split("\t", $HASH_repository_kindex{$kindex_K});
			$HASH_info_Kx{"k-mer"}			= $ARRAY_repository[1];
			$FILE_repository_info 			= $ARRAY_repository[4]."repository_".$kindex_K.".info";				
				
			#get all infos for repository
			my $href_info_Kx 	= &read_config($FILE_repository_info, \%HASH_info_Kx, \%HASH_repository_kindex, "run");
			$ARRAY_HASH_info_aref[$ki] 		= $href_info_Kx;
				
			if(!defined $compare){
				#START RUN SK
				%HASH_info_Kx				= %{$href_info_Kx};		
				&run_kmasker_SK($fasta, $kindex_K, \%HASH_info_Kx, \%HASH_repository_kindex);
			}					
		}
	
		if(defined $compare){
			#COMPARE
			&run_kmasker_MK($fasta, \@multi_kindex, \@ARRAY_HASH_info_aref, \%HASH_repository_kindex);		
		}
	}
	
	#QUIT
	print "\n - Thanks for using Kmasker! -\n\n";
		if($verbose) {
        	print "Output of external commands was written to " . $kmasker::kmasker_run::log."\n";
        }
    else{
        	unlink($kmasker::kmasker_run::log);
        }
    &createREADME();
	exit();
}


#######################
###
### EXPLORE SECTION
###
if(defined $explore){
	#USE EXPLORE MODULE
	
	########
	#STORE EXPLORE INFO in HASH info
	$HASH_info{"user name"}				= $user_name;
	$HASH_info{"kindex"}				= $kindex			if(defined $kindex);
	$HASH_info{"multi kindex"}			= \@multi_kindex 	if(defined $multi_kindex[0]);
	#ADDITIONAL
	$HASH_info{"version KMASKER"}		= $version;
	$HASH_info{"version BUILD"} 		= "";	
	
	print "\n starting explore module ... session id: $PID\n";
	if($verbose){
		if(defined $force) {
			print "\n --force is active!\n";
			#$force = "true";
		}
	}
	
	#ANNOTATION
	if(defined $custom_annotate){
		# EXPLORE ANNOTATE by GFF
		
		#START annotation of sequence features, provided in GFF with custome FASTA sequence or blastableDB
		my $check_settings = "";
		$check_settings .= " --gff undefined;" if(!defined $gff);
		$check_settings .= " --blast db or fasta undefined;" if((!defined $dbfasta)&&(!defined $blastableDB));
		$check_settings .= " --feature undefined;" if(!defined $feature);
		$check_settings .= " --fasta missing. It provides FASTA belonging to constructed GFF" if(!defined $fasta);
				
		if($check_settings eq ""){
    		$HASH_db{"db_fasta"} 	= $dbfasta if(defined $dbfasta);
    		$HASH_db{"db"} 			= $blastableDB if(defined $blastableDB);
			&custom_annotation($fasta, $gff, $feature ,\%HASH_db, \%HASH_info, \%HASH_path);
		}else{
			print "\n WARNING: Required parameter are missing.";
			print "\n WARNINGS are ".$check_settings;
			print "\n Kmasker was stopped.\n\n";
			exit;
		}		
	}	
	
	#REPORT
	if(defined $cstats){
		# EXPLORE STATS		
		#STATS requirements
		my $check_settings = 1;
		$check_settings = 0 if(!defined $occ);
		
		#CALL comparative methods
		if(scalar(@OCClist)>1){
			my $k_correct = 21;
			$k_correct = $k if(defined $k);	
			system("$path/OCC_compare.pl --k ".$k_correct." --fc ".$foldchange." --occ1 ".$OCClist[0]." --occ2 ".$OCClist[1]);
			exit();
		}else{
			print "\n You need to provide two OCC files. Process stopped!\n\n";
			exit();
		}
	}
	
	
	#REPORT
	if(defined $stats){
		# EXPLORE STATS
		
		#STATS requirements
		my $check_settings = 1;
		$check_settings = 0 if(!defined $occ);
		
		#exit
		if(scalar(@OCClist)>1){
			print "\n Too many files provided in '--occ'. Statistics only can be calculated for a single OCC file.\n\n";
			exit();
		}
				
		#CALL REPORT method
		$check_settings = 0 if(!defined $gff);		
		if($check_settings == 1){
    		&report_statistics($occ, $gff, $out);
    		
    		if(defined $max){
    			system("OCC_average.pl --occ ".$occ." --max ".$max." --out KMASKER_average_m".$max."_".$occ);
    		}else{
    			system("OCC_average.pl --occ ".$occ." --out KMASKER_average_".$occ);
    		}
		}else{
			print "\n WARNING: Missing required parameter '--occ' or '--gff'. Kmasker was stopped.\n\n";
			exit;
		}		
	}	
	
	
	#EXTRACT
	if(defined $xtract){
		# EXPLORE EXTRACT non-masked regions
		
		#STATS requirements
		my $check_settings = 1;
		$check_settings = 0 if(!defined $fasta);
				
		if($check_settings == 1){
    		kmasker::functions::Xtract($fasta, $length_threshold);
		}else{
			print "\n WARNING: Required parameter --fasta is missing. Kmasker was stopped.\n\n";
			exit;
		}		
	}
	
	#VIZUALISATION
	my $visualisation = 0;
	#simple plots
	$visualisation = 1 if(defined $hist);
	$visualisation = 1 if(defined $histm);
	$visualisation = 1 if(defined $barplot);
	#comparative plots
	$visualisation = 2 if(defined $violin);
	$visualisation = 2 if(defined $hexplot);
	$visualisation = 2 if(defined $boxplot);
		
	#HISTOGRAM
	if($visualisation == 1){
		# EXPLORE VIZUALISATIONS
		
		if(defined $occ){
		# explore visualisations require an OCC file
		
			my $missing_parameter = "";
			if(! -e $occ) {
				print "\n ERROR: $occ was not found. Kmasker was stopped.\n\n";
				exit;
			}
		
			#HISTOGRAM RAW or MEAN
			if((defined $hist)||(defined $histm)){	
				if(defined $list){				
					if (-e $list) {
   						&plot_histogram_raw($occ, $list, $force) if(defined $hist);
   						&plot_histogram_mean($occ, $list, $dynamic, $force, $sws) if(defined $histm);
					}else{
						$missing_parameter .= " --list (file ".$list." not found)";
					}				
				}else{
					#check
					my $systemcall_grep 		= "grep -c \">\" ".$occ;
					my $number 	= `$systemcall_grep`;
					$number		=~ s/\n//;
					if($number >= 10000){	
						
						if(defined $force){
							&plot_histogram_raw($occ, undef, $force) if(defined $hist);
							&plot_histogram_mean($occ, undef, $dynamic, $force, $sws) if(defined $histm);
						}else{						
							print "\n\n WARNING: Your input contains ".$number." sequences. We do not recommend to build more than 10.000 in one step.\n";
							print     " This might take very long. You can force this by using '--force'. Kmasker has been stopped \n\n";	
							$missing_parameter .= " --list";						
							exit();
						}
					}else{
							&plot_histogram_raw($occ, undef, $force) if(defined $hist);
							&plot_histogram_mean($occ, undef, $dynamic, $force, $sws) if(defined $histm);
					}
					#clean logs
					#if((! defined $verbose) && -e "log.txt") {
					#	unlink("log.txt");
					#}
				}
			}
			
			#BARPLOT
			if(defined $barplot){				
				my $outname = "KMASKER_average_".$occ;
				if(defined $max){
    				system("OCC_average.pl --occ ".$occ." --max ".$max." --out KMASKER_average_m".$max."_".$occ);
    				$outname = "KMASKER_average_m".$max."_".$occ; 
    			}else{
    				system("OCC_average.pl --occ ".$occ." --out ".$outname);
    			}	
    			$outname =~ s/occ$/txt/;			
				&plot_barplot($outname, $bin) if(defined $bin);
				&plot_barplot($outname) if(!defined $bin);
			}			
			
			if($missing_parameter ne ""){
				#GIVE warning note for missing parameter
				print "\n ERROR: missing parameter (".$missing_parameter.") !\n\n";
			}
			
		}else{
			print "\n ERROR: no occ provided. For this 'Kmasker --explore' subfunction an occ file is required!\n\n";
		}
	}
	# VISUALISATION 2 
	elsif($visualisation == 2){
		#COMPARATIVE PLOTS
		
		my $plottype = "";
		$plottype = "boxplot" if(defined $boxplot);
		$plottype = "violin"  if(defined $violin);
		$plottype = "hexplot" if(defined $hexplot);
		
		if(scalar(@OCClist)>1){
			#CALL comparative methods			
			#system("$path/OCC_compare.pl --fc ".$foldchange." --occ1 ".$OCClist[0]." --occ2 ".$OCClist[1]." --out KMASKER_comparative_descriptives_FC".$foldchange.".stats");
			#$cfile = "KMASKER_comparative_descriptives_FC".$foldchange.".stats";
			
			print "\n\ .. Please use '--cstats' before calling this method and subsequently provide result with '--cfile'. Process stopped!\n\n";
			exit();
		}
		
		my %HASH_list = ();
		if(defined $list){				
			if (-e $list) {
				#READING
				my $LIST 	 = new IO::File($list, "r") or die "\n unable to read list $list $!";
				while(<$LIST>){
					next if($_ =~ /^$/);
					next if($_ =~ /^#/);
					my $line = $_;
					$line =~ s/\n//;
					$line =~ s/ /\t/g;
					my @ARRAY_in = split(/\t/, $line);			
					$HASH_list{$ARRAY_in[0]} = 1;
				}   				
   			}else{
				print "\n .. cannot find file ".$list."\n";
				print "    process stopped!\n";
				exit();
   			}
		}		
		
		#FILE SIMPLE & DIRECT
		if(defined $file){
			if(-e $file){
				print "\n .. working with ".$file."\n";
			}else{
				print "\n .. cannot find file ".$file."\n";
				print "    process stopped!\n";
				exit();
			}
			#START	
			&plot_maker_direct($file, $plottype) if(defined $file);
		}
	
		
		#CFILE USER & COMPLEX
		if(defined $cfile){
			if(-e $cfile){
				print "\n .. working with ".$cfile."\n";
			}else{
				print "\n .. cannot find file ".$cfile."\n";
				print "    process stopped!\n";
				exit();
			}
			#START	
			&plot_maker($cfile, $plottype, \%HASH_list);
		}
		
	}
	
	#QUIT
	print "\n\n - Thanks for using Kmasker! -\n\n";
	if($verbose) {
        	print "Output of external commands was written to " . $kmasker::kmasker_explore::log ."\n";
        }
    else{
        	unlink($kmasker::kmasker_explore::log);
        }
    &createREADME();
	exit();
}
	

##END MAIN


## subroutine
#
sub show_repository(){	
	print "\n\nREPOSITORY of available kindex structures:\n";
	
	my %HASH_repository_content = ();
			
	#PRINT
	foreach my $kindex_this (keys %HASH_repository_kindex){
		my @ARRAY_repository	= split("\t", $HASH_repository_kindex{$kindex_this});
		my $hashline = "";
		
		$hashline = "\n\t";
		$hashline .= sprintf("%-25s", $kindex_this);				
		$hashline .= "\t".$ARRAY_repository[1]."\t";
		$hashline .= sprintf("%-14s", $ARRAY_repository[2]);
		$hashline .= "\t".$ARRAY_repository[3];
		
		#INSERT
		$HASH_repository_content{$kindex_this} = $hashline;
	}
	
	foreach my $key (sort keys %HASH_repository_content){
		print $HASH_repository_content{$key};
	}
	
	print "\n\n";
}

## subroutine
#
sub show_details_for_kindex(){
	my $kindex_this = $_[0];
	if(exists $HASH_repository_kindex{$kindex_this}){
		my @ARRAY_repository	= split("\t", $HASH_repository_kindex{$kindex_this});
		my $BUILD_file 		= new IO::File($ARRAY_repository[4]."repository_".$kindex_this.".info") or die " ... can not read/find repository_".$kindex_this.".info file in ".$ARRAY_repository[4]." details : $!\n\n";
		print "\n  KINDEX details for ".$kindex_this.": \n";
		while(<$BUILD_file>){
			print "\t".$_;
		}
	}else{
		print "\n\n WARNING: Requested kindex (".$kindex_this."). does not exist. Please check and use different index name.";
	}	
	print "\n\n";
}


## subroutine
#
sub initiate_user(){
	
	#$user_name 			= `whoami`;
	#$user_name			=~ s/\n//g;
	#my $uconf 			= $ENV{"HOME"}."/.kmasker_user.config";

	if(-e $uconf){
		#USER already exists, do nothing
		print "\n Using $uconf as user configuration. \n" if(defined $verbose);
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
			if(scalar(@seq_usr) == 0){
				print "\n .. kmasker was stopped: no input sequence provided (--seq) !";
				print "\n\n";
				exit(0);
			}
		}elsif(defined $set_external_path){
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
		
		if(((!defined $fasta)&&(!defined $krispr))){
			print "\n .. kmasker was stopped: no sequence provided (either '--fasta' or '--grna') !";
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

	#$user_name 			= `whoami`;
	#$user_name			=~ s/\n//g;
	#my $gconf 			= $path."/kmasker.config";
	
	#PERMISSION - calling this procedure is only be possible for directory owner (who installed Kmasker)
	#my $fp			 =  $path."/kmasker.config";
	print "\n Global configuration: $gconf\n" if(defined $verbose);
	if (-e $gconf) {
		my ($installed_by) = getpwuid(stat($gconf)->uid);
		#$installed_by =~ s/\n//;
		if($installed_by ne $user_name){
			print "\n Your user rights are not sufficient to create or modify this global configuration file.\n";
			print "I=(".$installed_by.") U=(".$user_name.")\n";
			exit();	
		}
	}
	
	#REQUIREMENTs
	my %HASH_requirments	= (	"jellyfish" => "",
								"fastq-stats" => "",
								"gffread" => "",
								"blastn"  => "",
								"makeblastdb" => "",
								"R" => "");
			
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
		my $gCFG_old 	= new IO::File($gconf, "r") or die "\n unable to read global config $!";	
		my $gCFG 		= new IO::File($gconf.".tmp", "w") or die "\n unable to update global config $!";
		
		my %HASH_provided = ();
		while(<$gCFG_old>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_tmp = split("=", $line);
			
			$HASH_provided{"jellyfish"} 	= $line if($line =~ /^jellyfish=/);
			$HASH_provided{"fastq-stats"} 	= $line if($line =~ /^fastq-stats=/);
			$HASH_provided{"gffread"}		= $line if($line =~ /^gffread=/);
			$HASH_provided{"blastn"} 	= $line if($line =~ /^blastn=/);
			$HASH_provided{"makeblastdb"} 	= $line if($line =~ /^makeblastdb=/);
			$HASH_provided{"R"}		= $line if($line =~ /^R=/);				
		}
		
		
		#CHECK tool requirments
		#JELLYFISH
		$HASH_requirments{"jellyfish"} = &check_routine_for_requirement("jellyfish", $HASH_provided{"jellyfish"}, $HASH_requirments{"jellyfish"});
		
		#FASTQ-STATs
		$HASH_requirments{"fastq-stats"} = &check_routine_for_requirement("fastq-stats", $HASH_provided{"fastq-stats"}, $HASH_requirments{"fastq-stats"});
		
		#GFFREAD
		$HASH_requirments{"gffread"} = &check_routine_for_requirement("gffread", $HASH_provided{"gffread"}, $HASH_requirments{"gffread"});

		#blastn
		$HASH_requirments{"blastn"} = &check_routine_for_requirement("blastn", $HASH_provided{"blastn"}, $HASH_requirments{"blastn"});
		
		#makeblastdb
		$HASH_requirments{"makeblastdb"} = &check_routine_for_requirement("makeblastdb", $HASH_provided{"makeblastdb"}, $HASH_requirments{"makeblastdb"});
		
		#R
		$HASH_requirments{"R"} = &check_routine_for_requirement("R", $HASH_provided{"R"}, $HASH_requirments{"R"});
				
		#WRITE
		print $gCFG "#external tool requirements\n";
		foreach my $required (keys %HASH_requirments){
			if($required !~ /^PATH_kindex/){
				system("which $HASH_requirments{$required} >/dev/null 2>&1 || { echo >&2 \"Kmasker requires $required but it's not installed or path is missing! Kmasker process stopped.\"; exit 1; \}");
				print $gCFG $required."=".$HASH_requirments{$required}."\n";
				print "\n info ".$required." --> ".$HASH_requirments{$required};
			}			
		}
		
		print "\n\n";
		$gCFG_old->close();
		$gCFG->close();
		#system("mv ".$gconf.".tmp ".$gconf)	
		move "$gconf.tmp", $gconf;

	} else { #global config is missing
		my $gCFG 		= new IO::File($gconf, "w") or die "\n unable to create global config $!";
		#CHECK tool requirments
		#JELLYFISH
		$HASH_requirments{"jellyfish"} = &check_routine_for_requirement("jellyfish", "", $HASH_requirments{"jellyfish"});
		
		#FASTQ-STATs
		$HASH_requirments{"fastq-stats"} = &check_routine_for_requirement("fastq-stats", "", $HASH_requirments{"fastq-stats"});
		
		#GFFREAD
		$HASH_requirments{"gffread"} = &check_routine_for_requirement("gffread", "", $HASH_requirments{"gffread"});

		#blastn
		$HASH_requirments{"blastn"} = &check_routine_for_requirement("blastn", "", $HASH_requirments{"blastn"});
		
		#makeblastdb
		$HASH_requirments{"makeblastdb"} = &check_routine_for_requirement("makeblastdb", "", $HASH_requirments{"makeblastdb"});
		
		#R
		$HASH_requirments{"R"} = &check_routine_for_requirement("R", "", $HASH_requirments{"R"});
				
		#WRITE
		print $gCFG "#external tool requirements\n";
		foreach my $required (keys %HASH_requirments){
			if($required !~ /^PATH_kindex/){
				system("which $HASH_requirments{$required} >/dev/null 2>&1 || { echo >&2 \"Kmasker requires $required but it's not installed or path is missing! Kmasker process stopped.\"; exit 1; \}");
				print $gCFG $required."=".$HASH_requirments{$required}."\n";
				print "\n info ".$required." --> ".$HASH_requirments{$required};
			}			
		}
		
		print "\n\n";
		$gCFG->close();
	}
}


## subroutine
#
sub use_expert_settings(){
	
	#INPUT has type: par1=value1;par2=value2; ... 
	
	#HASH_info: par1 value; par2 value; ...
	
	my $type 		= $_[0];
	my $parameter 	= $_[1];
	
	#PARSING
	my @ARRAY_user_parameter_input = split(";", $parameter);
	
	if($type eq "blast"){		
		print "\n\n BLAST expert settings should be provided in the form '-parameter:value;-parameter:value'";
		print "\n e.g. blast parameter '-perc_identity 0.95 -evalue 10' encoded for Kmasker '-perc_identity=0.95;-evalue=10' ";
		print "\n Please note that the usage of expert settings will overwrite ALL blast default settings (set those manually, too)!";
		print "\n For reusability you can use the '--config_blast' parameter to provide a file with your manual settings for blast.\n\n";
	}	
	
	my $hash_info_entry = "";
	
	foreach(@ARRAY_user_parameter_input){
		if($_ =~ /=/){
			my @ARRAY_tmp = split("=", $_);
			my $PAR = $ARRAY_tmp[0];
			if(defined $ARRAY_tmp[1]){
				my $VALUE = $ARRAY_tmp[1];
				
				#KMASKER
				if($type eq "kmasker"){
					$hash_info_entry = "";
					$hash_info_entry = $HASH_info{"user setting kmasker"} if(exists $HASH_info{"user setting kmasker"}); 
					
					if($PAR eq "pctgap"){
						print "\n\n .. changing default parameter for PCTGAP from ".$MK_percent_gapsize." to ". $VALUE." !\n";
						$MK_percent_gapsize = $VALUE;
						$HASH_info{"MK_percent_gapsize"}	= $MK_percent_gapsize;
						$hash_info_entry					.= " " if($hash_info_entry ne "");
						$hash_info_entry 					.= "--pctgap ".$MK_percent_gapsize.";";
					}
					
					if($PAR eq "minseed"){
						print "\n\n .. changing default parameter for MINSEED from ".$MK_min_seed." to ". $VALUE." !\n";
						$MK_min_seed = $VALUE;
						$HASH_info{"MK_min_seed"}			= $MK_min_seed;
						$hash_info_entry					.= " " if($hash_info_entry ne "");
						$hash_info_entry 					.= "--minseed ".$MK_min_seed.";";
					}
					
					if($PAR eq "mingff"){
						print "\n\n .. changing default parameter for MINGFF from ".$MK_min_gff." to ". $VALUE." !\n";
						$MK_min_gff 						= $VALUE;
						$HASH_info{"MK_min_gff"}			= $MK_min_gff;
						$hash_info_entry					.= " " if($hash_info_entry ne "");
						$hash_info_entry 					.= "--mingff ".$MK_min_gff.";";
					}
					
					if($PAR eq "rept"){
						print "\n\n .. changing default parameter for REPT from ".$repeat_threshold." to ". $VALUE." !\n";
						$repeat_threshold					= $VALUE;
						$HASH_info{"rept"}					= $repeat_threshold;
						$hash_info_entry					.= " " if($hash_info_entry ne "");
						$hash_info_entry 					.= "--rept ".$repeat_threshold.";";
					}
					
					if($PAR eq "minl"){
						print "\n\n .. changing default parameter for MINL from ".$length_threshold." to ". $VALUE." !\n";
						$length_threshold 					= $VALUE;
						$HASH_info{"minl"}					= $length_threshold;
						$hash_info_entry					.= " " if($hash_info_entry ne "");
						$hash_info_entry 					.= "--minl ".$length_threshold.";";
					}
					
					if($PAR eq "bed"){
						print "\n\n .. changing default parameter for BED from ".$HASH_info{"bed"}." to ". $VALUE." !\n";
						$HASH_info{"bed"}					= $VALUE;
						$hash_info_entry					.= " " if($hash_info_entry ne "");
						$hash_info_entry 					.= "--bed ".$VALUE.";";
					}
										
					#INSERT
					$HASH_info{"user setting kmasker"} 		= $hash_info_entry;
				}	
				
				
				#JELLYFISH
				if($type eq "jelly"){
					$hash_info_entry = "";
					$hash_info_entry = $HASH_info{"user setting jelly"} if(exists $HASH_info{"user setting jelly"}); 
					if($PAR eq "size"){
						print "\n .. changing default parameter for --size from ".$size." to ". $VALUE." !\n";
						$size 				= $VALUE;
						$hash_info_entry	.= " " if($hash_info_entry ne "");
						$hash_info_entry 	.= "--size ".$size;
					}
					
					if($PAR eq "threads"){
						print "\n .. changing default parameter for --threads from ".$threads." to ". $VALUE." !\n";
						$threads			= $VALUE;
						$hash_info_entry	.= " " if($hash_info_entry ne "");
						$hash_info_entry 	.= "--threads ".$threads;
					}
					
					#INSERT
					$HASH_info{"user setting jelly"} 		= $hash_info_entry;
				}
				
				
				#BLAST
				if($type eq "blast"){					
					## blast parameter encoded for Kmasker '-perc_identity=0.95;-evalue=10' ";
					$hash_info_entry = "";
					$hash_info_entry = $HASH_info{"user setting blast"} if(exists $HASH_info{"user setting blast"}); 
					$PAR = "-".$PAR if($PAR !~ /^\-/);					
					$hash_info_entry	.= " " if($hash_info_entry ne "");
					$hash_info_entry 	.= $PAR." ".$VALUE;
			
					#INSERT
					$HASH_info{"user setting blast"} 		= $hash_info_entry;					
				}			
				
			}else{
				print "\n\n Cannot use user input (".$parameter."). Kmasker is stopped!\n\n\n";
				exit();
			}
		}else{
			print "\n\n Cannot use user input (".$parameter."). Kmasker is stopped!\n\n\n";
			exit();
		}
	}
	
	print "\n .. modified parameter settings for ".$type." : ".$hash_info_entry."\n" if($hash_info_entry ne "");
	
}


## subroutine
#
sub read_configuration_file(){	
	
	my $type = $_[0];
	my $file = $_[1];
	
	if(-e $file){
	
		if($type eq "blast"){		
			print "\n\n READING expert settings from configuration file ".$file." for blast.\n"; 	
			my $parsed = &reading_file($type, $file);
			print "\n USER SETTING blast: '".$parsed."'\n";
			$HASH_info{"user setting blast"} 		= $parsed;
		}
				
		if($type eq "jelly"){		
			print "\n\n READING expert settings from configuration file ".$file." for jellyfish.\n"; 
			my $parsed = &reading_file($type, $file);
			print "\n USER SETTING jelly: '".$parsed."'\n";
			$HASH_info{"user setting jelly"} 		= $parsed;
		}
				
		if($type eq "kmasker"){		
			print "\n\n READING expert settings from configuration file ".$file." for kmasker.\n"; 
			my $parsed = &reading_file($type, $file);
			print "\n USER SETTING Kmasker: '".$parsed."'\n";
			$HASH_info{"user setting kmasker"} 		= $parsed;
		}
	}else{
		print "\n Configuration file '".$file."' is not found!\n";
	}	
}


sub reading_file(){
	
	my $type = $_[0];
	my $file = $_[1];
		
	my $CONFIG = new IO::File($file, "r") or die "\n unable to read user config ".$file." $!";	
	my $readline = "";	
	my @ARRAY_collect = ();
	while(<$CONFIG>){
		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		my $line = $_;
		$line =~ s/\n//;
		$line =~ s/^ //g;
		push(@ARRAY_collect, $line);
	}
	$CONFIG->close();	
	$readline 			= join(" ", @ARRAY_collect);
	my $readline_func 	= join(";", @ARRAY_collect);
	$readline_func 		= join("=", @ARRAY_collect) if($type eq "blast");
	&use_expert_settings($type, $readline_func);
		
	#RETURN	
	return $readline;
}

sub check_tool () {
my $line = $_[0];
my $tool = $_[1]; 
if($line =~ /^$tool=/){
	my @ARRAY_tmp = split("=", $line);
	if(!defined $ARRAY_tmp[1]){
		system("which $tool >/dev/null 2>&1 || { echo >&2 \"Kmasker requires $tool but it's not installed! Kmasker process stopped.\"; exit 1; \}");
		$HASH_path{"$tool"} = `which $tool`;
		$HASH_path{"$tool"} =~ s/\n//;
		}else{
			$HASH_path{"$tool"} = $ARRAY_tmp[1];
			system("which ".$HASH_path{"$tool"}." >/dev/null 2>&1 || { echo >&2 \"Kmasker requires $tool but it's not installed!  Kmasker process stopped.\"; exit 1; \}");
		}
		print "\n $tool=".$HASH_path{"$tool"}."\n" if(defined $verbose);
	}
}

## subroutine
#
sub read_user_config(){
	#$user_name 			= `whoami`;
	#$user_name			=~ s/\n//g;
	#my $gconf 			= $path."/kmasker.config";
	#my $uconf 			= $ENV{"HOME"}."/.kmasker_user.config";
	
	if(-e $gconf){
		#LOAD info for external tools
		my $gCFG = new IO::File($gconf, "r") or die "\n unable to read global config $!";	
		
		while(<$gCFG>){
			next if($_ =~ /^$/);
			next if($_ =~ /^#/);
			my $line = $_;
			$line =~ s/\n//;
			my @ARRAY_tmp = split("=", $line);
			if($ARRAY_tmp[0] eq "PATH_kindex_external"){
				if(defined $ARRAY_tmp[1]){
					$PATH_kindex_external	= $ARRAY_tmp[1];
					$PATH_kindex_external .= "/" if($PATH_kindex_external !~ /\/$/);
					
					#EXTERNAL
					if(-d $PATH_kindex_external){
						#directory exists - do nothing
					}else{
						#directory has to be created
						system("mkdir "."\"".$PATH_kindex_external."\"");
					}
				}
			}		
			
			#READ external tool path
			&check_tool($line, "jellyfish");
			&check_tool($line, "fastq-stats");
			&check_tool($line, "gffread");
			&check_tool($line, "blastn");
			&check_tool($line, "makeblastdb");
			&check_tool($line, "R");

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

			&check_tool($line, "jellyfish");
			&check_tool($line, "fastq-stats");
			&check_tool($line, "gffread");
			&check_tool($line, "blastn");
			&check_tool($line, "makeblastdb");
			&check_tool($line, "R");
			
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
			
			#EXTERNAL			FIXME ... does not fix the problem
			if(($ARRAY_tmp[0] =~ "export")&&(defined $ARRAY_tmp[1])){
				system($ARRAY_tmp[0]."=".$ARRAY_tmp[1]);
				if(defined $verbose){
					print "\n user settings applied:\n";
					print " ".$ARRAY_tmp[0]."=".$ARRAY_tmp[1]."\n";
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
	if($PATH_kindex_external !~ /^(\/)$|^(\s)$|^()$/){
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
	
}

#################
#
#
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
				my $path_which = `which $requirement`;
				$path_which =~ s/\n$//;
				if($path_which eq "") {
					print "\n No success! Your configuration for $requirement is might be wrong!\n";
					print "\n Current path is $default!\n";
				} else {
					print "\n Success! Set path for $requirement to $path_which!\n";
					$default = $path_which;
			    }

			}else{
				$default = $path_check;
			}					
		}
	}	
	return $default;
}


#################
#
#
#ub create_PID(){
	#FIXME
	# @Chris: please continue to create PID and LONG_PID
	
#	my $loctime = localtime;
#	$loctime = strftime('%Y%m%d%H%M%S',localtime); ## outputs 120817100834	
	
#	$PID		= $loctime;	#PLACE HOLDER
#	$LONG_PID	= $loctime;	
#}


#################
#
#
sub intro_call(){

	#SHOW version
	if(defined $show_version){
		print "\n Kmasker version: ".$version."\n\n";
		exit();
	}

	#CHECK settings
	if(defined $check_install){
		&check_install();
		exit();
	}

	#CHECK export
	my $export			= $path."/.kmasker_settings.sh";
	if(-e $export){
		system(".".$export);
		print "\n EXPORT paths\n\n";
	}
	
	#CHECK existence of FASTA
	if(defined $fasta){
		if(!(-e $fasta)){
			print "\n INPUT FASTA \'".$fasta."\' does not exists!";	
			print "\n .. process stopped! \n\n";		
			exit();	
		}
		my $cwd = getcwd();
		my $filename = fileparse($fasta);
		if(!(-e $cwd . "/" . $filename)) {
			symlink($fasta, $cwd . "/" . $filename);
			$fasta = $filename;
		}
		$fasta=fileparse(fasta_to_uppercase($fasta));


	}
	
	#CHECK existence of SEQ
	if((scalar @seq_usr) > 0){
		my $check_status = 1;
		foreach my $seq (@seq_usr){
			if(!(-e $seq)){
				print "\n INPUT seq \'".$seq."\' does not exists!";		
				$check_status = 0;	
			}
		}
		if($check_status == 0){
			print "\n .. process stopped! \n\n";		
			exit();
		}
	}
	
	#SET private path
	if(defined $set_private_path){
		
		if($set_private_path =~ /^\./){
			print "\n .. please use absolut path! Kmasker was stopped!\n";
			exit();
		}
		
		print "\n Kmasker will change private path: ".$set_private_path;
		print "\n Please note that data will not be moved to new directories!";
		print "\n This has to be done manually!\n\n";
		
		&set_private_path($set_private_path);
		exit();
	}
	
	#SET external path
	if(defined $set_external_path){
		
		if($set_external_path =~ /^\./){
			print "\n .. please use absolut path! Kmasker was stopped!\n";
			exit();
		}
		
		print "\n Kmasker will change external path: ".$set_external_path;
		print "\n Please note that data will not be moved to new directories!";
		print "\n This has to be done manually!\n\n";
		
		&set_external_path($set_external_path);
		exit();
	}
	
	#GET info of Kmasker used path
	if(defined $show_path){
		&show_path_infos();
		exit();
	}	
	
	#GENERAL options
	#SHOW repository
	if(defined $show_kindex_repository){
		&show_repository();
		exit();
	}	
		
	#SHOW details	
	if(defined $show_details_for_kindex){
		
		&show_details_for_kindex($show_details_for_kindex);
		exit();
	}
	
	#REMOVE kindex
	if(defined $remove_kindex){
		&remove_kindex($remove_kindex,\%HASH_repository_kindex);
		exit();
	}
	
	#OCClist
	if(scalar(@OCClist) == 1){
		$occ = $OCClist[0];
	}	
	
	#THREADS
	if(defined $foldchange_usr){
		$foldchange = $foldchange_usr;
	}
	
	#MISMATCH
	if(defined $mismatch_usr){
		$mismatch = $mismatch_usr;
	}

}


#################
#
#
sub help(){

	print "\n Usage of program Kmasker:";
    print "\n (version:  ".$version.")\t\t\t\t(session id: $PID)";
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
		print "\n --as \t\t input is of sequence type 'assembly'. Set this option if your input is based on ";
		print "\n      \t\t assembled sequences (e.g. genome reference). Default sequence type is 'reads'.";
		
		print "\n\n";
		exit();
	}
	
	if(defined $run){
		#HELP section run
		print "\n Command:";
		print "\n\t Kmasker --run --fasta sequence_to_be_analyzed.fasta\n";		
		
		print "\n\n Options:";
		print "\n --kindex\t single or multiple k-mer indices (use space delimited list)";
		print "\n --fasta\t sequences in FASTA format for k-mer analysis and masking";
		print "\n --krispr\t provide gRNA sequences in FASTA format for specificity analysis";
		print "\n --compare\t perform comparative analysis using multiple k-mer indices (requires --kindex K1 K2)";
		print "\n --rept\t\t frequency threshold used for masking [5]!";
		print "\n --minl\t\t minimal length of sequence. Kmasker will extract all non-repetitive sequences with sufficient length [100]";
	
		print "\n\n";
		exit();
	}
	
	if(defined $explore){
		#HELP section explore
		print "\n Command (subset):";
		print "\n\n\t Kmasker --explore --annotate --fasta query.fasta --gff kmasker_result.gff --feature KRC --dbfasta repeats.fasta";
		print "\n\n\t Kmasker --explore --hist --occ file.occ --list sequence.ids";
		print "\n\n\t Kmasker --explore --hexplot --kindex At1 Hv1";
		print "\n\n\t Kmasker --explore --stats --occ file.occ";
		
		print "\n\n Options:";
		print "\n --annotate\t\t custom annotation using featured elements of GFF (requires --gff, --fasta, --db or --dbfasta)";
		print "\n --gff\t\t\t use Kmasker constructed GFF report for annotation";
		print "\n --feature\t\t the type of feature in the GFF that should be annotated";
		print "\n --dbfasta\t\t custom sequences [FASTA] with annotated features in sequence descriptions";
		print "\n --db\t\t\t pre-calculated blastableDB of nucleotides used for annotation";	
		
		print "\n --xtract\t\t extract non-masked regions from Xmasked FASTA (requires --fasta)";
		print "\n --hist\t\t\t create histogram using raw values (requires --occ and optional --list)";
		print "\n --histm\t\t create histogram using calulated means (requires --occ and optional --list)";
		print "\n --violin\t\t create violin plot (apply after '--run --compare' mode)";
		print "\n --hexplot\t\t create hexagon plot (apply after '--run --compare' mode)";
		print "\n --barplot\t\t create barplot (apply after '--run --compare' mode)";
		print "\n --occ\t\t\t provide a Kmasker constructed occ file containing k-mer frequencies";
		print "\n --file\t\t\t provide simple file (three columns) for plotting";
		print "\n --cfile\t\t\t provide comparative file for plotting";
		print "\n --list\t\t\t file containing a subset of selected contig identifier";	
		print "\n --stats\t\t create report of basic statistics (requires --occ and --gff)";	
		print "\n --cstats\t\t create copmarative statistics (requires --occ)";
		
		
		print "\n\n";
		exit();
	}
	

    print "\n Description:\n\t Kmasker is a tool for the automatic detection of repetitive sequence regions.";
    print "\n\t There are three modules and you should select one for your analysis.";
    
    print "\n\n Modules:";
	print "\n --build\t\t construction of new index (requires --indexfiles)";
	print "\n --run\t\t\t perform analysis and masking (requires --fasta)";
	print "\n --explore\t\t perform downstream analysis with constructed index and detected repeats";
	
	print "\n\n General options:";
	print "\n --show_repository\t\t show complete list of private and external k-mer indices";
	print "\n --show_details\t\t\t show details for a requested kindex";
	print "\n --show_path\t\t\t show path Kmaskers looks for constructed kindex";
	print "\n --remove_kindex\t\t remove kindex from repository";
	print "\n --set_private_path\t\t change path to private repository";
	print "\n --set_external_path\t\t change path to external repository [readonly]";
	print "\n --expert_setting_kmasker\t submit individual parameter to Kmasker eg. pctgap,";
	print "\n\t\t\t\t minseed, mingff (see documentation!)";
	print "\n --expert_setting_jelly\t\t submit individual parameter to jellyfish (e.g. on memory usage ";
	print "\n\t\t\t\t for index construction)";
	print "\n --expert_setting_blast\t\t submit individual parameter to blast (e.g. '-evalue')";
	print "\n --threads\t\t\t set number of threads [4]";
	print "\n --bed\t\t\t\t force additional BED output [off]";
	print "\n --user_conf\t\t set specific user configuration file [$uconf]";
	print "\n --global_conf\t\t set specific global configuration file [$gconf]";
	print "\n --long_pid\t\t create a process id that is unique for this host (e.g. for use in cluster environments)";


	
	print "\n\n";
	exit();
}

sub createREADME {
   open(my $readme, ">", "KMASKER_00_README_$PID.txt") or die "Can not open KMASKER_00_README_$PID.txt\n";
   print $readme "PID: $PID\n";
   print $readme "Kmasker version: $version\n";
   print $readme "Input command line: ";
   print $readme join(" ", @ARGV_B);
   print $readme "\n";
   print $readme "Mode: ";
   if(defined $run){
   		print $readme "RUN ";
   		if(defined $krispr) {
   			print $readme "KRISPR ";
   		}
   		if(defined $fish){
   			print $readme "FISH ";
   		}
   		if(scalar(@multi_kindex)==1){
   			print $readme "REPEAT_MASKING ";
   		}
   		if(scalar(@multi_kindex)>1){
   			print $readme "COMPARE";
   		}
   }
   if(defined $explore) {
   		print $readme "EXPLORE ";
   		if(defined $custom_annotate){
   			print $readme "ANNOTATE ";
   		}
   		if(defined $cstats){
   			print $readme "COMPARE_STATS ";
   		}
   		if(defined $stats){
   			print $readme "STATS ";
   		}
   		if(defined $xtract){
   			print $readme "EXTRACT "
   		}
   		if((defined $barplot) || (defined $hist) || (defined $histm) || (defined $violin)){
   			print $readme "VISUALISATION ";
   		}
   		if(defined $hist){
   			print $readme "HISTOGRAM_RAW ";
   		}
   		if(defined $histm){
   			print $readme "HISTOGRAM_MEAN "; 
   		}
   }

   print $readme "\nInput files: ";
   foreach(($fasta, $krispr, $gff, $dbfasta ,$occ, $list)) {
   	if(defined $_) {
   		print $readme "\n\t\t" . $_;
   	}
   }
   print $readme "\n";
   print $readme "Output files ";
   
    my $dir = '.';

    opendir(DIR, $dir) or die $!;

    my @files 
        = grep { 
            /$PID/             
	} readdir(DIR);

    # Loop through the array printing out the filenames
    foreach my $file (@files) {
        print $readme "\n\t\t$file";
    }

    closedir(DIR);
    
    foreach(("kmasker_plots_$PID", "kmasker_raw_plots_$PID", "$temp_path")) {
	    if(-e $_) {  
		    my $dir = $_;

		    opendir(DIR, $dir);

		    my @temp = grep { 
            /($PID)|(png$)/             
			} readdir(DIR);

		    # Loop through the array printing out the filenames
		    foreach my $file (@temp) {
		        print $readme "\n\t\t$dir/$file";
		    }
		    closedir(DIR);
		}
	}
    print $readme "\n";	
}



1;
