package kmasker::kmasker_run;
use Exporter qw(import);
use File::Basename;
use strict;
use warnings;
use kmasker::filehandler;
use kmasker::occ;
use kmasker::functions;
use File::Basename qw(dirname);
use Cwd  qw(abs_path);

#adapt
our @ISA = qw(Exporter);
our @EXPORT = qw(
	run_kmasker_SK
	run_kmasker_MK
	show_version_PM_run
);
our @EXPORT_OK = qw(run_kmasker_SK run_kmasker_MK run_gRNA show_version_PM_run);


## VERSION
my $version_PM_run 	= "0.0.30 rc180802";

## subroutine
#
sub run_kmasker_SK{
	# SINGLE KINDEX (SK)
	my $fasta 			= $_[0];
	my $kindex			= $_[1];
	my $href_info		= $_[2];
	my $href_repo 		= $_[3];
	my $path 		= dirname abs_path $0;		
	
	my %HASH_info_this 			= %{$href_info};
	my %HASH_repository_kindex 	= %{$href_repo};
   	
	if(exists $HASH_repository_kindex{$kindex}){
		
		# GET info from repository
		my $rept 				= $HASH_info_this{"rept"};
		my $length_threshold 	= $HASH_info_this{"min_length"};
		my $seq_depth			= sprintf("%.0f" , $HASH_info_this{"sequencing depth"});
		my $k					= $HASH_info_this{"k-mer"};
		my $threads				= $HASH_info_this{"threads"};
		my $temp_path       	= $HASH_info_this{"temp_path"};
		my $verbose				= $HASH_info_this{"verbose"};
		my $percent 			= $HASH_info_this{"MK_percent_gapsize"}; 	# default : 10 (%)
		my $min_seed			= $HASH_info_this{"MK_min_seed"};			# default :  5 (bp)
		my $min_gff				= $HASH_info_this{"MK_min_gff"}; 			# default : 10 (bp)
		my $bed					= $HASH_info_this{"bed"}; 					# default : 'no'
			
		my @ARRAY_repository	= split("\t", $HASH_repository_kindex{$kindex});
		my $absolut_path		= $ARRAY_repository[4];

		print "\n parameter setting: rept       = ".$rept;
		print "\n parameter setting: min_length = ".$length_threshold;
		print "\n parameter setting: minseed    = ".$min_seed;
		print "\n parameter setting: pctgap     = ".$percent;
		print "\n parameter setting: mingff     = ".$min_gff;
		print "\n parameter setting: bed     	= ".$bed;
		print "\n temp path: $temp_path";
		print "\n\n";

		#create symbolic link to kindex from private or global
		my $full_kindex_name = "KINDEX_".$kindex.".jf";		
		if(-e $absolut_path.$full_kindex_name){
			system("ln -s \"".$absolut_path.$full_kindex_name."\"");
		}else{
			print "\n WARNING: KINDEX (".$full_kindex_name.") not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
			exit();
		}
		
		#start
		system("$path/cmasker -f \"".$fasta."\" -j \"".$full_kindex_name."\" -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex . " >>log.txt 2>&1");
       
       #clean
		unlink($full_kindex_name);
        if(!(-e "KMASKER_".$kindex."_RT".$rept."_NORM"."_".$fasta)) {
            print "\n .. KMASKER_".$kindex."_RT".$rept."_NORM"."_".$fasta." was not generated!\n";
            print "\n .. please provide a bug report!\n\n";
            exit();
        }
        mkdir($temp_path, 0775);
       
        #make tab from masked fasta
        print "\n .. start to generate TAB" ;#if(!defined $silent);
        system ("mv" . " *.occ \"$temp_path\"");
        kmasker::filehandler::fasta_to_tab("KMASKER_".$kindex."_RT".$rept."_NORM"."_".$fasta, $temp_path); #just change .fasta to .tab in temp
        kmasker::filehandler::sequence_length($fasta);
        system( "mv" ." $fasta.length" . " \"$temp_path/$fasta.length\"" );
        
		#merge seeds
		print "\n .. start to generate extended region" ;#if(!defined $silent);
		my $tab = $fasta;
		$tab =~ s/(\.fasta$)|(\.fa$)//; 
		kmasker::filehandler::merge_tab_seeds("$temp_path/KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab.".tab", $percent, $min_seed);
		
		#PRODUCE GFF
		print "\n .. start to generate GFF" ;#if(!defined $silent);
		my $feature = "KRC";
		my $subfeature = "KRR";
		kmasker::filehandler::tab_to_gff("$temp_path/KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab."_Regions_merged.tab", "$temp_path/$fasta.length" ,$min_gff, $feature ,"$temp_path/KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab.".tab", $subfeature);
		
		#extract non-repetitive regions
		kmasker::functions::Xtract("KMASKER_".$kindex."_RT".$rept."_NORM"."_".$fasta, $length_threshold);
       
       	my $gffname = "RESULT_KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab.".gff";
        system("mv" . " $temp_path/KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab."_Regions_merged.gff " . $gffname);
        
        if($bed eq "1"){
        	#Feedback
			print "\n .. start to generate BED" ;#if(!defined $silent);
        	#BED format output is activated
        	&write_gff2bed($gffname); 		
        }
        
        #Statistics
        #Feedback
        print "\n .. start to generate statistics" ;#if(!defined $silent);
        system("$path/stats.R " . "-i " . "\"" .$temp_path . "\"" . "/KMASKER_" . $kindex . "_NORM_" . $tab . ".occ" . " -g " . "RESULT_KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab.".gff" . " -c sequence" . " >>log.txt 2>&1");
        system("$path/stats.R "  . "-i " . "\"" .$temp_path . "\"" . "/KMASKER_" . $kindex . "_NORM_" . $tab . ".occ" . " -g " . "RESULT_KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab.".gff" . " -c " . $feature . " >>log.txt 2>&1");
        if((!(-e "Report_statistics_". $feature . "_RESULT_KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab.".gff.tab")) || (!(-e  "Report_statistics_sequence_" . "RESULT_KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab. ".gff.tab" ))) {
        	print "Some statistics could not be calculated. The main reason for this is that there are no significant features in the gff file.\n";
        }
        else {
	       	system("$path/stats_overview.R " . " -s Report_statistics_sequence_" . "RESULT_KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab. ".gff.tab " . " -k " . "Report_statistics_". $feature . "_RESULT_KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab.".gff.tab" . " >>log.txt 2>&1");
        	system("mv " . "overview_stats.txt " . "Report_overview.txt")
        }
        if($verbose) {
        	print "Output of external commands was written to log.txt\n";
        }
        else{
        	unlink("log.txt");
        }

	}else{
		#KINDEX is missing in repository
		print "\n .. Kmasker was stopped!\n";
		print "\n .. The kindex (".$kindex.") you requested is not available in repository!\n\n";
	}	
}


## subroutine
#  run kmasker with multiple kindex (MK) structures for comparative analysis
sub run_kmasker_MK{
	# MULTIPLE KINDEX (MK)
	my $fasta 			= $_[0];
	my $aref_mkindex	= $_[1];
	my $aref_info_Kx	= $_[2];
	my $href_repo 		= $_[3];
	#my $href_info		= $_[4];
	my $path 			= dirname abs_path $0;
	my $FC_compare		= 5;
		
	
	#my %HASH_info_this 			= %{$href_info};
	my %HASH_repository_kindex 	= %{$href_repo};
	
	my @ARRAY_kindex 			= @{$aref_mkindex};
	
	# GET info from repository
	my @ARRAY_INFO_Kx		= @{$aref_info_Kx};
	my %HASH_info_this		= %{$ARRAY_INFO_Kx[0]};
		
	# GET info
	my $rept 				= $HASH_info_this{"rept"};
	my $length_threshold 	= $HASH_info_this{"min_length"};
	my $percent 			= $HASH_info_this{"MK_percent_gapsize"}; 	# default : 10 (%)
	my $min_seed			= $HASH_info_this{"MK_min_seed"};			# default :  5 (bp)
	my $min_gff				= $HASH_info_this{"MK_min_gff"}; 			# default : 10 (bp)
	my $bed					= $HASH_info_this{"bed"}; 					# default : 'no'
	my @ARRAY_full_kindex_names = ();
	my @ARRAY_seq_depth			= ();
	my $temp_path       	= $HASH_info_this{"temp_path"};
	my $verbose				= $HASH_info_this{"verbose"};
	
	print "\n parameter setting: rept       = ".$rept;
	print "\n parameter setting: min_length = ".$length_threshold;
	print "\n parameter setting: minseed    = ".$min_seed;
	print "\n parameter setting: pctgap     = ".$percent;
	print "\n parameter setting: mingff     = ".$min_gff;
	print "\n parameter setting: bed     	= ".$bed;
	print "\n temp path: $temp_path";
	print "\n";
	
	#ADD info
	my @ARRAY_occ = ();
	
	#USER info
	my $fold_change = 10;
	
	for(my $i=0;$i<scalar(@ARRAY_kindex);$i++){
		
		my $kindex 				= $ARRAY_kindex[$i];   		
		my $aref_kindexhash		= $ARRAY_INFO_Kx[$i];
   		my %HASH_info_this		= %{$aref_kindexhash};
   		my $seq_depth			= $HASH_info_this{"sequencing depth"};
		my $k					= $HASH_info_this{"k-mer"};	
		my $threads				= $HASH_info_this{"threads"};
		print "\n .. starting multi kindex processing on ".$kindex."\n\n";
   		   		
   		#READ repository.info   		
		if(exists $HASH_repository_kindex{$kindex}){		
			
			my @ARRAY_repository 	= split("\t", $HASH_repository_kindex{$kindex});
			my $absolut_path		= $ARRAY_repository[4];

			#create symbolic link to kindex from private or global
			my $full_kindex_name = "KINDEX_".$kindex.".jf";
			push(@ARRAY_full_kindex_names, $full_kindex_name);
			push(@ARRAY_seq_depth, $seq_depth);
			if(-e $absolut_path.$full_kindex_name){
				system("ln -s ".$absolut_path.$full_kindex_name);
				
				#PRODUCE OCC
				system("$path/cmasker -f ".$fasta." -j ".$full_kindex_name." -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex." -s" . " >>log.txt 2>&1");
				
				my @NAME = split(/\./, $fasta);
				pop @NAME;
				my $name = join(".", @NAME);
				print "\n NAME = ".$name."\n";
				push(@ARRAY_occ, "KMASKER_".$kindex."_NORM"."_".$name.".occ");
				
				#clean
    			system("rm ".$full_kindex_name);
				
			}else{
				print "\n ... using path ".$absolut_path;
				print "\n WARNING: KINDEX (".$full_kindex_name.") not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
				exit();
			}
		}else{
			#KINDEX is missing in repository
			print "\n .. Kmasker was stopped!\n";
			print "\n .. The kindex (".$kindex.") you requested is not available in repository!\n\n";
		}	
	} ## close foreach LOOP	
	
    mkdir($temp_path, 0775);
	kmasker::filehandler::sequence_length($fasta);
    system( "mv" ." $fasta.length" . " \"$temp_path/$fasta.length\"" );
	#Feedback
	print "\n .. start to generate TAB" ;#if(!defined $silent);
		
	#Convertion to TAB file
	my $occ1 = $ARRAY_occ[0];
	my $occ2 = $ARRAY_occ[1];
	kmasker::occ::multi_occ($rept, $fold_change, $occ1, $occ2, "$temp_path/KMASKER_comparativ_FC".$fold_change."_");
	(my $name1,my $path1,my $suffix1) = fileparse($occ1, qr/\.[^.]*/);
	(my $name2,my $path2,my $suffix2) = fileparse($occ2, qr/\.[^.]*/);
	
	#Modify here to support more than two files!
	#Feedback
	print "\n .. start to generate extended region" ;#if(!defined $silent);
	
	#EXTEND regions (for TAB1 and TAB2)	
	my $tab1 = $occ1;
	my $tab2 = $occ2;
	$tab1 =~ s/\.occ$//;
	$tab2 =~ s/\.occ$//;
	kmasker::filehandler::merge_tab_seeds("$temp_path/KMASKER_comparativ_FC".$fold_change."_".$tab1.".tab", $percent, $min_seed);
	kmasker::filehandler::merge_tab_seeds("$temp_path/KMASKER_comparativ_FC".$fold_change."_".$tab2.".tab", $percent, $min_seed);
	
	#Feedback
	print "\n .. start to generate GFF" ;#if(!defined $silent);
	
	#PRODUCE GFF
	my $feature 	= "KRC";
	my $subfeature 	= "KRR";
	my $gffname_D1 	= "KMASKER_comparativ_FC".$fold_change."_$tab1" . ".gff";
	my $gffname_D2 	= "KMASKER_comparativ_FC".$fold_change."_$tab2" . ".gff";	
	kmasker::filehandler::tab_to_gff("$temp_path/KMASKER_comparativ_FC".$fold_change."_$tab1" . "_Regions_merged.tab", "$temp_path/$fasta.length", $min_gff, $feature, "$temp_path/KMASKER_comparativ_FC".$fold_change."_".$tab1.".tab", $subfeature);
	kmasker::filehandler::tab_to_gff("$temp_path/KMASKER_comparativ_FC".$fold_change."_$tab2" . "_Regions_merged.tab", "$temp_path/$fasta.length", $min_gff, $feature, "$temp_path/KMASKER_comparativ_FC".$fold_change."_".$tab2.".tab", $subfeature);
	system("cp" . " $temp_path/KMASKER_comparativ_FC".$fold_change."_$tab1" . "_Regions_merged".".gff " . $gffname_D1);
    system("cp" . " $temp_path/KMASKER_comparativ_FC".$fold_change."_$tab2" . "_Regions_merged".".gff " . $gffname_D2);
    
    if($bed eq "1"){
    	#Feedback
		print "\n .. start to generate BED" ;#if(!defined $silent);
   		#BED format output is activated
    	&write_gff2bed($gffname_D1); 
    	&write_gff2bed($gffname_D2); 		
    }
    
    #Statistics
    system("$path/stats.R " . "-i " . $occ1 . " -g"  .$gffname_D1. " -c sequence" .  " >>log.txt 2>&1");
    system("$path/stats.R " . "-i " . $occ2 . " -g " .$gffname_D2. " -c sequence" .  " >>log.txt 2>&1");

    system("$path/stats.R " . "-i " . $occ1 . " -g"  .$gffname_D1. " -c " . $feature . " >>log.txt 2>&1");
    system("$path/stats.R " . "-i " . $occ2 . " -g " .$gffname_D2. " -c " . $feature . " >>log.txt 2>&1");

    
    if( (!(-e "Report_statistics_sequence_" . "KMASKER_comparativ_FC".$fold_change."_$tab1" . ".gff.tab")) || (!(-e "Report_statistics_". $feature . "_KMASKER_comparativ_FC".$fold_change."_$tab1" . ".gff.tab"))) {
        	print "\nSome statistics could not be calculated (" .  $ARRAY_kindex[0] ."). The main reason for this is that there are no significant features in the gff file.\n";
        }
    else {
 		system("$path/stats_overview.R " . " -s Report_statistics_sequence_" . "KMASKER_comparativ_FC".$fold_change."_$tab1" . ".gff.tab" . " -k " . "Report_statistics_". $feature . "_KMASKER_comparativ_FC".$fold_change."_$tab1" . ".gff.tab " . " >>log.txt 2>&1");
 		system("mv overview_stats.txt " . "Report_overview_" . $ARRAY_kindex[0] . ".txt");
     }
   # system("$path/stats_overview.R " . " -s Stats_Sequence_" . "KMASKER_comparativ_FC".$fold_change."_$tab1" . ".gff.tab" . " -k " . "Stats_". $feature . "KMASKER_comparativ_FC".$fold_change."_$tab1" . ".gff.tab " . " >>log.txt 2>&1");
    if(!(-e "Report_statistics_sequence_" . "KMASKER_comparativ_FC".$fold_change."_$tab2" . ".gff.tab") || !(-e "Report_statistics_". $feature . "_KMASKER_comparativ_FC".$fold_change."_$tab2" . ".gff.tab")) {
        	print "\nSome statistics could not be calculated (" .  $ARRAY_kindex[1] ."). The main reason for this is that there are no significant features in the gff file.\n";
        }
    else {
 		system("$path/stats_overview.R " . " -s Report_statistics_sequence_" . "KMASKER_comparativ_FC".$fold_change."_$tab2" . ".gff.tab" . " -k " . "Report_statistics_". $feature . "_KMASKER_comparativ_FC".$fold_change."_$tab2" . ".gff.tab " . " >>log.txt 2>&1");
 		system("mv overview_stats.txt " . "Report_overview_" . $ARRAY_kindex[0] . ".txt");
     }

    #CALL comparative methods
    system("$path/OCC_compare.pl --fc ".$FC_compare." --occ1 ".$occ1." --occ2 ".$occ2."");
    if($verbose) {
      	print "Output of external commands was written to log.txt\n";
    }
    else{
       	unlink("log.txt");
    }
}

## subroutine
#  run gRNA will check characteristics of custome FASTA
sub run_gRNA(){
	my $kindex_this = $_[0];
	my $gRNA		= $_[1];
	my $href_repo	= $_[2];
		
	#GET INFO
	my $path 					= dirname abs_path $0;	
	my %HASH_repository_kindex 	= %{$href_repo};
	my @ARRAY_repository 		= split("\t", $HASH_repository_kindex{$kindex_this});
	my $absolut_path			= $ARRAY_repository[4];
		
	print "\n\n ... start Kmasker gRNA module\n";	
	my $full_kindex_name = "KINDEX_".$kindex_this.".jf";		
	if(-e $absolut_path.$full_kindex_name){
		system("ln -s \"".$absolut_path.$full_kindex_name."\"");
	}else{
		print "\n WARNING: KINDEX (".$full_kindex_name.") not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
		exit();
	}	
	system("$path/crispr.py ".$kindex_this." ".$gRNA);	
	print "\n\n ... Kmasker gRNA module finished \n";
}


## subroutine
#
sub get_kindex_info{
	my $href_info_this = $_[0];
	my $href_repo_this = $_[1];
	my %HASH_info_tmp = %{$href_info_this};
	my %HASH_repo_tmp = %{$href_repo_this};
	### not required
	
}

## subroutine
#
sub show_version_PM_run{
	print "\n\nVERSION of module run: ".$version_PM_run."\n\n";
}

## subroutine
#
sub write_gff2bed{
	my $gffname = $_[0];
	
	my $bedname = $gffname;
    $bedname =~ s/\.gff$/.bed/;
    my $GFFFILE = new IO::File($gffname, "r") or die "\n unable to read gff ".$gffname." $!";	
    my $BEDFILE = new IO::File($bedname, "w") or die "\n unable to write bed ".$bedname." $!";
        	
    #WRITE
    while(<$GFFFILE>){
    	next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		my @ARRAY_gff = split("\t", $_);
		my $substring = "KRC";
		if($ARRAY_gff[2] =~ m/$substring/){
			print $BEDFILE $ARRAY_gff[0]."\t".$ARRAY_gff[3]."\t".$ARRAY_gff[4]."\n";
		}
	}
		   	
	#CLOSE
	$GFFFILE->close();
	$BEDFILE->close();	
}

	
	

1;
