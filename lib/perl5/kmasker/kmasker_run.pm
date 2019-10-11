package kmasker::kmasker_run;
use Exporter qw(import);
use File::Basename;
use File::Copy;
use strict;
use warnings;
use kmasker::filehandler;
use kmasker::occ;
use kmasker::functions;
use File::Basename qw(dirname);
use Cwd  qw(getcwd abs_path);

my $timestamp = getLoggingTime();
our $log = "log_run_" . $kmasker::functions::PID . ".txt";

#adapt
our @ISA = qw(Exporter);
our @EXPORT = qw(
	run_kmasker_SK
	run_kmasker_MK
	show_version_PM_run
	$log
);
our @EXPORT_OK = qw(run_kmasker_SK run_kmasker_MK run_krispr show_version_PM_run create_krispr_model $log);


## VERSION
my $version_PM_run 	= "0.0.35 rc190212";


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
		my $length_threshold 	= $HASH_info_this{"minl"};
		my $seq_depth			= sprintf("%.0f" , $HASH_info_this{"sequencing depth"});
		my $k					= $HASH_info_this{"k-mer"};
		my $threads				= $HASH_info_this{"threads"};
		my $temp_path       	= $HASH_info_this{"temp_path"};
		my $verbose				= $HASH_info_this{"verbose"};
		my $percent 			= $HASH_info_this{"MK_percent_gapsize"}; 	# default : 10 (%)
		my $min_seed			= $HASH_info_this{"MK_min_seed"};			# default :  5 (bp)
		my $min_gff				= $HASH_info_this{"MK_min_gff"}; 			# default : 10 (bp)
		my $bed					= $HASH_info_this{"bed"}; 					# default : 'no'
		my $PID					= $HASH_info_this{"PID"};
        my $strict              = $HASH_info_this{"strict"};
			
		my @ARRAY_repository	= split("\t", $HASH_repository_kindex{$kindex});
		my $absolut_path		= $ARRAY_repository[4];
		print "\n .. starting kindex processing on ".$kindex."\n\n";
		print "\n parameter setting: rept       = ".$rept;
		print "\n parameter setting: minl       = ".$length_threshold;
		print "\n parameter setting: minseed    = ".$min_seed;
		print "\n parameter setting: pctgap     = ".$percent;
		print "\n parameter setting: mingff     = ".$min_gff;
		print "\n parameter setting: bed        = ".$bed;
		print "\n temp path: $temp_path";
		print "\n\n";

		#create symbolic link to kindex from private or global
		my $full_kindex_name = "KINDEX_".$kindex.".jf";	
		if(-e $absolut_path.$full_kindex_name){
			my $sl_result = eval {symlink("${absolut_path}${full_kindex_name}", getcwd()."/".$full_kindex_name); 1};
					if (($sl_result == 0) ||! (-e $full_kindex_name)) {
						die "Symbolic link of ${absolut_path}${full_kindex_name} could not be created!\n";
					}	
		}else{
			print "\n WARNING: KINDEX (".$full_kindex_name.") not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
			exit();
		}
		
		# OCC NAME
		my $occ_kmer_counts = "KMASKER_kmer_counts_KDX_".$kindex."_".$PID.".occ";		
		my $masked_fasta;
        my $strict_param = "";
        if(defined $strict){
            $strict_param = " -t ";
        }
		#start
		system("cmasker -f \"".$fasta."\" -j \"".$full_kindex_name."\" -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex . "_" . $PID . $strict_param . " >>$log 2>&1");
       
       #clean
		unlink($full_kindex_name);
        if(!(-e "KMASKER_masked_KDX_".$kindex."_".$PID.".fasta")) {
            print "\n .. KMASKER_masked_KDX_".$kindex."_".$PID.".fasta" . " was not generated!\n";
            print "\n .. please provide a bug report!\n\n";
            exit();
        }
        else {
        	$masked_fasta = "KMASKER_masked_KDX_".$kindex."_".$PID.".fasta";
        }
        
        #BED #FIXME
		if($bed eq "1"){
		   	#BED format output is activated
		  	my @HELP = split(".", $fasta);
		   	pop @HELP;
		   	my $name = join(".", @HELP);		   	
		   	print "\n NAME = ".$name;		   	
		   	system("OCC2BED.pl --occ $occ_kmer_counts >>$log 2>&1");		
   		}
   		        
        mkdir($temp_path, 0775);
       
        #make tab from masked fasta
        print "\n .. start to generate TAB file" ;#if(!defined $silent);
        #system ("mv" . " *.occ \"$temp_path\"");
        kmasker::filehandler::fasta_to_tab($masked_fasta, $temp_path); #just change .fasta to .tab in temp
        kmasker::filehandler::sequence_length($masked_fasta);
        move("${masked_fasta}.length" , "$temp_path/${masked_fasta}.length" );
        
		#MERGE SEEDS
		print "\n .. start to generate extended region" ;#if(!defined $silent);
		my $tab = $masked_fasta;
		$tab =~ s/(\.fasta$)|(\.fa$)//; 
		kmasker::filehandler::merge_tab_seeds("$temp_path/$tab.tab", $percent, $min_seed);
		
		#PRODUCE GFF
		print "\n .. start to generate GFF" ;#if(!defined $silent);
		my $feature = "KRC";
		my $subfeature = "KRR";
		kmasker::filehandler::tab_to_gff("$temp_path/$tab"."_regions_merged.tab", "$temp_path/${masked_fasta}.length" ,$min_gff, $feature ,"$temp_path/$tab". ".tab", $subfeature);
		
		#extract non-repetitive regions
		kmasker::functions::Xtract($masked_fasta, $length_threshold, "KDX_${kindex}_$PID");
       
       	my $gffname = "KMASKER_repeat_regions_KDX_${kindex}_$PID.gff";
        move("$temp_path/${tab}_regions_merged.gff",$gffname);
        
        if($bed eq "1"){
        	#Feedback
			print "\n .. start to generate BED" ;#if(!defined $silent);
        	#BED format output is activated
        	kmasker::filehandler::write_gff2bed($gffname); 		
        }
        
        #Statistics
        #Feedback
        print "\n .. start to generate statistics\n" ;#if(!defined $silent);
        system("$path/stats.R " . "--input " . $occ_kmer_counts . " --gff " . $gffname . " --class sequence" . " --pid $PID" . " >>$log 2>&1");
        system("$path/stats.R " . "--input " . $occ_kmer_counts . " --gff " . $gffname . " --class " . $feature . " --pid $PID" . " >>$log 2>&1");
        if((!(-e "KMASKER_report_statistics_". $feature . "_$PID.tab")) || (!(-e  "KMASKER_report_statistics_sequence" . "_$PID.tab"))) {
        	print "\n Some statistics could not be calculated. The main reason for this is that there are no significant features in the gff file.\n";
        }
        else {
        	move("KMASKER_report_statistics_". $feature . "_$PID.tab", "KMASKER_report_statistics_". $feature . "_$PID" . ".tab");
        	move("KMASKER_report_statistics_". "sequence" . "_$PID.tab", "KMASKER_report_statistics_". "sequence" . "_$PID" . ".tab");
	       	system("$path/stats_overview.R " . " -s KMASKER_report_statistics_". "sequence" . "_$PID" . ".tab" . " -k " . "KMASKER_report_statistics_". $feature . "_$PID" . ".tab" . " -p $PID" . " >>$log 2>&1");
	       	if(!-e("KMASKER_report_overview_statistics_$PID.txt")) {
	       		print "\n Statistics overview could not be calculated!\n";
	       	}
	       	else {
	       		move("KMASKER_report_overview_statistics_$PID.txt", "KMASKER_report_overview_$PID.txt");
	       	}
        }
	}else{
		#KINDEX is missing in repository
		print "\n .. Kmasker was stopped!\n";
		print "\n .. The kindex (".$kindex.") you requested is not available in the repositories!\n\n";
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
	my $length_threshold 	= $HASH_info_this{"minl"};
	my $percent 			= $HASH_info_this{"MK_percent_gapsize"}; 	# default : 10 (%)
	my $min_seed			= $HASH_info_this{"MK_min_seed"};			# default :  5 (bp)
	my $min_gff				= $HASH_info_this{"MK_min_gff"}; 			# default : 10 (bp)
	my $bed					= $HASH_info_this{"bed"}; 					# default : 'no'
	my @ARRAY_full_kindex_names = ();
	my @ARRAY_seq_depth			= ();
	my $temp_path       	= $HASH_info_this{"temp_path"};
	my $verbose				= $HASH_info_this{"verbose"};
	my $fold_change 		= $HASH_info_this{"fold-change"};
	my $PID					= $HASH_info_this{"PID"};
    my $strict              = $HASH_info_this{"strict"};


	print "\n parameter setting: rept       = ".$rept;
	print "\n parameter setting: minl       = ".$length_threshold;
	print "\n parameter setting: minseed    = ".$min_seed;
	print "\n parameter setting: pctgap     = ".$percent;
	print "\n parameter setting: mingff     = ".$min_gff;
	print "\n parameter setting: bed     	= ".$bed;
	print "\n temp path: $temp_path";
	print "\n";
	
	#ADD info
	my @ARRAY_occ 	= ();
	my $global_k;
	
	
	for(my $i=0;$i<scalar(@ARRAY_kindex);$i++){
		
		my $kindex 				= $ARRAY_kindex[$i];   		
		my $aref_kindexhash		= $ARRAY_INFO_Kx[$i];
   		my %HASH_info_this		= %{$aref_kindexhash};
   		my $seq_depth			= $HASH_info_this{"sequencing depth"};
		my $k					= $HASH_info_this{"k-mer"};	
		my $threads				= $HASH_info_this{"threads"};
		
		if(!defined $global_k){
			$global_k = $k;
		}else{
			if($global_k ne $k){
				print "\n WARNING: You are using/comparing KINDEX with differnt k-mer size!!! This is not recommended!!!\n\n";
			}
		}
		
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
				my $sl_result = eval {symlink("${absolut_path}${full_kindex_name}", getcwd()."/".$full_kindex_name); 1};
				if (($sl_result == 0) ||! (-e $full_kindex_name)) {
					die "Symbolic link of ${absolut_path}${full_kindex_name} could not be created!\n";
				}
				
				#PRODUCE OCC
                my $strict_param = "";
                 if(defined $strict){
                     $strict_param = " -t ";
                 }
                #start
                system("cmasker -f \"".$fasta."\" -j \"".$full_kindex_name."\" -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex . "_" . $PID . $strict_param . " >>$log 2>&1");
				#my @NAME = split(/\./, $fasta);
				#pop @NAME;
				#my $name = join(".", @NAME);
				#print "\n NAME = ".$name."\n";
				push(@ARRAY_occ, "KMASKER_kmer_counts_KDX_".$kindex."_$PID.occ");
				
				#BED
				 if($bed eq "1"){ #FIX ME
			    	#BED format output is activated
			    	system("OCC2BED.pl --occ KMASKER_kmer_counts_KDX_".$kindex."_$PID.occ --rept ".$rept." >>$log 2>&1");		
   				}
				
				#clean
    			unlink($full_kindex_name);
				
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
    move("$fasta.length" , "$temp_path/$fasta.length");
	#Feedback
	print "\n .. start to generate TAB" ;#if(!defined $silent);
		
	#Convertion to TAB file
	my $occ1 = $ARRAY_occ[0];
	my $occ2 = $ARRAY_occ[1];
	kmasker::occ::multi_occ($rept, $fold_change, $occ1, $occ2, "$temp_path/KMASKER_comparative_FC".$fold_change."_${PID}_");
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
	kmasker::filehandler::merge_tab_seeds("$temp_path/KMASKER_comparative_FC".$fold_change."_${PID}_".$tab1.".tab", $percent, $min_seed);
	kmasker::filehandler::merge_tab_seeds("$temp_path/KMASKER_comparative_FC".$fold_change."_${PID}_".$tab2.".tab", $percent, $min_seed);
	
	#Feedback
	print "\n .. start to generate GFF" ;#if(!defined $silent);
	
	#PRODUCE GFF
	my $feature 	= "KDC";
	my $subfeature 	= "KDR";
	my $gffname_D1 	= "KMASKER_diverse_regions_KDX_" .$ARRAY_kindex[0].  "_$PID" . ".gff";
	my $gffname_D2 	= "KMASKER_diverse_regions_KDX_" .$ARRAY_kindex[1].  "_$PID" . ".gff";
	kmasker::filehandler::tab_to_gff("$temp_path/KMASKER_comparative_FC".$fold_change."_${PID}_$tab1" . "_regions_merged.tab", "$temp_path/$fasta.length", $min_gff, $feature, "$temp_path/KMASKER_comparative_FC".$fold_change."_${PID}_".$tab1.".tab", $subfeature);
	kmasker::filehandler::tab_to_gff("$temp_path/KMASKER_comparative_FC".$fold_change."_${PID}_$tab2" . "_regions_merged.tab", "$temp_path/$fasta.length", $min_gff, $feature, "$temp_path/KMASKER_comparative_FC".$fold_change."_${PID}_".$tab2.".tab", $subfeature);
	move("${temp_path}/KMASKER_comparative_FC".$fold_change."_${PID}_$tab1" . "_regions_merged".".gff" , getcwd() . "/". $gffname_D1);
    move("${temp_path}/KMASKER_comparative_FC".$fold_change."_${PID}_$tab2" . "_regions_merged".".gff" , getcwd() . "/". $gffname_D2);
    
    if($bed eq "1"){ #FIXME 
    	#Feedback
		print "\n .. start to generate BED" ;#if(!defined $silent);
   		#BED format output is activated
    	kmasker::filehandler::write_gff2bed($gffname_D1, $feature); 
    	kmasker::filehandler::write_gff2bed($gffname_D2, $feature); 		
    }
    
    #Statistics
    print "\n .. start to generate statistics\n" ;#if(!defined $silent);

    system("$path/stats.R " . "--input " . $occ1 . " --gff "  .$gffname_D1. " --class sequence" . " --pid $PID" . " >>$log 2>&1");
    system("$path/stats.R " . "--input " . $occ1 . " --gff "  .$gffname_D1. " --class " . $feature . " --pid $PID" ." >>$log 2>&1");

    if((!(-e "KMASKER_report_statistics_sequence_$PID.tab")) || (!(-e "KMASKER_report_statistics_${feature}_$PID.tab"))) {  	
      	print "\nSome statistics could not be calculated (" .  $ARRAY_kindex[0] ."). The main reason for this is that there are no significant features in the gff file.\n";
    }
    else {
 		system("$path/stats_overview.R " . " -s KMASKER_report_statistics_sequence_$PID.tab" . " -k " . "KMASKER_report_statistics_". $feature . "_$PID.tab " . " -p $PID" . " >>$log 2>&1");
 		move("KMASKER_report_statistics_sequence_$PID.tab", "KMASKER_report_statistics_KDX_". $ARRAY_kindex[0] ."_sequence_$PID.tab");
 		move("report_overview_statistics_$PID.txt" , "KINDEX_report_overview_statistics_KDX_" . $ARRAY_kindex[0] . "_$PID.txt");
 	}

 	system("$path/stats.R " . "--input " . $occ2 . " --gff "  .$gffname_D2. " --class sequence" . " --pid $PID" . " >>$log 2>&1");
    system("$path/stats.R " . "--input " . $occ2 . " --gff "  .$gffname_D2. " --class " . $feature . " --pid $PID" ." >>$log 2>&1");

    if((!(-e "KMASKER_report_statistics_sequence_$PID.tab")) || (!(-e "KMASKER_report_statistics_${feature}_$PID.tab"))) {  	
      	print "\nSome statistics could not be calculated (" .  $ARRAY_kindex[1] ."). The main reason for this is that there are no significant features in the gff file.\n";
    }
    else {
 		system("$path/stats_overview.R " . " -s KMASKER_report_statistics_sequence_$PID.tab" . " -k " . "KMASKER_report_statistics_". $feature . "_$PID.tab " . " -p $PID". " >>$log 2>&1");
 		move("KMASKER_report_statistics_sequence_$PID.tab", "KMASKER_report_statistics_KDX_". $ARRAY_kindex[1] ."_sequence_$PID.tab");
 		move("report_overview_statistics_$PID.txt" , "KINDEX_report_overview_statistics_KDX_" . $ARRAY_kindex[1] . "_$PID.txt");
 	}


    #CALL comparative methods
    system("$path/OCC_compare.pl --k ".$global_k." --fc ".$FC_compare." --occ1 ".$occ1." --occ2 ".$occ2." --out KMASKER_report_statistics_compare_$PID.txt");
}

## subroutine
#  run krispr will check characteristics of custome FASTA
sub run_krispr(){
	my $kindex_this = $_[0];
	my $krispr		= $_[1];
	my $href_repo	= $_[2];
	my $href_info	= $_[3];
	
	#GET INFO
	my $path 					= dirname abs_path $0;	
	my %HASH_repository_kindex 	= %{$href_repo};
	my %HASH_info_this		 	= %{$href_info};
	my @ARRAY_repository 		= split("\t", $HASH_repository_kindex{$kindex_this});
	my $absolut_path			= $ARRAY_repository[4];
	my $kripr_coverage_threshold= $HASH_info_this{"rept"};
	my $kripr_mismatch			= $HASH_info_this{"krisp mismatch"};
	my $PID						= $HASH_info_this{"PID"};
	my $model 					= $HASH_info_this{"krisp model"};
	#my $OUT_krispr				= "KMASKER_krispr_results.txt";
	my $OUT_krispr				= "KMASKER_krispr_KDX_".$kindex_this."_".$PID.".txt";
	
		
	#SETUP
	if(defined $model) {
		copy($model, ".") or die "Copy failed: $!";
	}
	else {
		copy($path."/../krispr/data_krispr.RData", "data_krispr.RData") or die "Copy failed: $!";
	}
	copy($path."/../krispr/models_krispr.R", ".") or die "Copy failed: $!";
	copy($path."/../krispr/krispr.py", ".") or die "Copy failed: $!";		
	print "\n\n ... start Kmasker krispr module\n";	
	my $full_kindex_name = "KINDEX_".$kindex_this.".jf";		
	if(-e $absolut_path.$full_kindex_name){
			my $sl_result = eval {symlink("${absolut_path}${full_kindex_name}", getcwd()."/".$full_kindex_name); 1};
					if (($sl_result == 0) ||! (-e $full_kindex_name)) {
						die "Symbolic link of ${absolut_path}${full_kindex_name} could not be created!\n";
					}	
	}else{
		print "\n WARNING: KINDEX (".$full_kindex_name.") not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
		exit();
	}
	
	#SINGLE SEQ
	# not activated in Kmasker
	# system("python3 ".$path."/krispr.py single -q ".$krispr_sequence." -j ".$full_kindex_name." -m ".$kripr_mismatch." -c ".$kripr_coverage_threshold);
	
	#MULTI FASTA
	system("python3 "."krispr.py multi -q ".$krispr." -j ".$full_kindex_name." -m ".$kripr_mismatch." -c ".$kripr_coverage_threshold." -t >".$OUT_krispr . ">>$log 2>&1");
	print "\n\n ... Kmasker krispr module finished \n";
	
	#CLEAN
	unlink($full_kindex_name);
	unlink("data_krispr.RData");
	unlink("models_krispr.R");
}
## subroutine
#  make new model for krispr
sub create_krispr_model(){
	my $path 					= dirname abs_path $0;	
	my $targets = $_[0];
	my $coverage		= $_[1]; 
	copy($path."/../krispr/data_krispr.RData", ".") or die "Copy failed: $!";
	copy($path."/../krispr/models_krispr.R", ".") or die "Copy failed: $!";
	copy($path."/../krispr/krispr.py", ".") or die "Copy failed: $!";
	copy($path."/../krispr/save_model_data.R", ".") or die "Copy failed: $!";
    system("python3 "."krispr.py new-model -e ".$targets." -c ".$coverage.">>$log 2>&1");
    unlink("data_krispr_backup.RData");
    if(-e "data_krispr.RData" ) {
    	print("\n\n A new model was created in the current directory. You can copy/rename it.\n You can use it in kmasker run with -model.\n");
    }
    unlink("models_krispr.R");
    unlink("krispr.py");
    unlink("save_model_data.R");
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

	
	

1;
