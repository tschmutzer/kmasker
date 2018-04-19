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
our @EXPORT_OK = qw(run_kmasker_SK run_kmasker_MK show_version_PM_run);


## VERSION
my $version_PM_run 	= "0.0.28 rc180419";

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
		my $percent 			= $HASH_info_this{"MK_percent_gapsize"}; 	# default : 10 (%)
		my $min_seed			= $HASH_info_this{"MK_min_seed"};			# default :  5 (bp)
		my $min_gff				= $HASH_info_this{"MK_min_gff"}; 			# default : 10 (bp)
			
		my @ARRAY_repository	= split("\t", $HASH_repository_kindex{$kindex});
		my $absolut_path		= $ARRAY_repository[4];

		print "\n parameter setting: rept       = ".$rept;
		print "\n parameter setting: min_length = ".$length_threshold;
		print "\n parameter setting: minseed    = ".$min_seed;
		print "\n parameter setting: pctgap     = ".$percent;
		print "\n parameter setting: mingff     = ".$min_gff;
		print "\n temp path: $temp_path";
		print "\n\n";

		#create symbolic link to kindex from private or global
		my $full_kindex_name = "KINDEX_".$kindex."_k".$k.".jf";		
		if(-e $absolut_path.$full_kindex_name){
			system("ln -s \"".$absolut_path.$full_kindex_name."\"");
		}else{
			print "\n WARNING: KINDEX (".$full_kindex_name.") not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
			exit();
		}
		
		#start
		system("$path/cmasker -f \"".$fasta."\" -j \"".$full_kindex_name."\" -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex);
       
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
		
		#Add annotation
		# TASK: Move to explore!
		#kmasker::functions::add_annotation($fasta, "$temp_path/KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab."_Regions_merged.tab", $BLASTDB, "$temp_path/KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab."_Regions_merged.gff", $threads);
        #system("FASTA_Xdivider.pl --fasta KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta." --sl ".$length_threshold);
        kmasker::functions::Xtract("KMASKER_".$kindex."_RT".$rept."_NORM"."_".$fasta);
        #system("mv Xsplit_KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta." KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_MIN".$length_threshold."_".$fasta);
        #system("mv" . " temp/Xsplit_KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta . " Xsplit_KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta);
        system("mv" . " $temp_path/KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab."_Regions_merged.gff" . " RESULT_KMASKER_".$kindex."_RT".$rept."_NORM"."_".$tab.".gff");
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
	my $path 		= dirname abs_path $0;		
	
	#my %HASH_info_this 			= %{$href_info};
	my %HASH_repository_kindex 	= %{$href_repo};
	
	my @ARRAY_kindex 			= @{$aref_mkindex};
	
	# GET info from repository
	my @ARRAY_INFO_Kx		= @{$aref_info_Kx};
	my %HASH_info_task		= %{$ARRAY_INFO_Kx[0]};
		
	# GET info
	my $rept 				= $HASH_info_task{"rept"};
	my $length_threshold 	= $HASH_info_task{"min_length"};
	my $percent 			= $HASH_info_this{"MK_percent_gapsize"}; 	# default : 10 (%)
	my $min_seed			= $HASH_info_this{"MK_min_seed"};			# default :  5 (bp)
	my $min_gff				= $HASH_info_this{"MK_min_gff"}; 			# default : 10 (bp)
	my @ARRAY_full_kindex_names = ();
	my @ARRAY_seq_depth			= ();
	my $temp_path       	= $HASH_info_task{"temp_path"};
	
	print "\n parameter setting: rept       = ".$rept;
	print "\n parameter setting: min_length = ".$length_threshold;
	print "\n parameter setting: minseed    = ".$min_seed;
	print "\n parameter setting: pctgap     = ".$percent;
	print "\n parameter setting: mingff     = ".$min_gff;
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

				
			# #REPEAT analytics
			# my $BLASTDB_redat 	= "";
			# my $BLASTDB_repbase = "";
			# my $BLASTDB_trep	= "";
			# my $BLASTDB = $repeat_lib_path ."/repeats.fasta";
			
			# #REPEAT analytics -should be done with REdat, RepBase and TREP FASTA
			# if(-e $repeat_lib_path."/REdat/mipsREdat_9.3p_ALL.fasta"){
			# 	$BLASTDB_redat		= "\"".$repeat_lib_path."/REdat/mipsREdat_9.3p_ALL.fasta\""; #CHANGEME
			# 	print "Using REdat: " . $BLASTDB_redat . "\n";
			# }
			
			# if(-e $repeat_lib_path."/RepBase/RepBase22.07_plants.fasta"){
			# 	my $BLASTDB_repbase	= "\"".$repeat_lib_path."/RepBase/RepBase22.07_plants.fasta\""; #CHANGEME
			# 	print "Using RepBase: " . $BLASTDB_repbase . "\n";

			# }
			
			# if(-e $repeat_lib_path ."/TREP/trep-db_nr_Rel-16.fasta"){
			# 	$BLASTDB_trep		= "\"".$repeat_lib_path."/TREP/trep-db_nr_Rel-16.fasta\""; #CHANGEME
			# 	print "Using TREP: " . $BLASTDB_trep . "\n";

			# }
			
			# if((! -e $repeat_lib_path."/repeats.fasta") || (-z $repeat_lib_path."/repeats.fasta")) {
			# 	system("cat ".$BLASTDB_redat." ".$BLASTDB_repbase." ".$BLASTDB_trep." > " .$repeat_lib_path  . "/repeats.fasta");
			# }
			# if(((! -e $repeat_lib_path."/repeats.fasta.nhr") || (-z $repeat_lib_path."/repeats.fasta.nhr")) || ((! -e $repeat_lib_path."/repeats.fasta.nin") || (-z $repeat_lib_path."/repeats.fasta.nin")) || ((! -e $repeat_lib_path."/repeats.fasta.nsq") || (-z $repeat_lib_path."/repeats.fasta.nsq"))) {
			# 	print("BLASTDB is missing...rebuilding it!\n");
			# 	system("makeblastdb -in \"".$BLASTDB."\" -dbtype nucl ");
			# 	print("BLASTDB was rebuilt.\n");
			# }
			# print "Using BLASTDB: " . $BLASTDB . "\n";


			#create symbolic link to kindex from private or global
			my $full_kindex_name = "KINDEX_".$kindex."_k".$k.".jf";
			push(@ARRAY_full_kindex_names, $full_kindex_name);
			push(@ARRAY_seq_depth, $seq_depth);
			if(-e $absolut_path.$full_kindex_name){
				system("ln -s ".$absolut_path.$full_kindex_name);
				
				#PRODUCE OCC
				system("cmasker -f ".$fasta." -j ".$full_kindex_name." -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex." -s");
				
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
	my $feature = "KRC";
	my $subfeature = "KRR";
	kmasker::filehandler::tab_to_gff("$temp_path/KMASKER_comparativ_FC10_$tab1" . "_Regions_merged.tab", "$temp_path/$fasta.length", $min_gff, $feature, "$temp_path/KMASKER_comparativ_FC10_".$tab1.".tab", $subfeature);
	kmasker::filehandler::tab_to_gff("$temp_path/KMASKER_comparativ_FC10_$tab2" . "_Regions_merged.tab", "$temp_path/$fasta.length", $min_gff, $feature, "$temp_path/KMASKER_comparativ_FC10_".$tab2.".tab", $subfeature);
	#Annotate GFF --> moved to explore
	#kmasker::functions::add_annotation($fasta, "temp/KMASKER_comparativ_FC10_$tab1" . "_Regions_merged.tab", $BLASTDB, "temp/KMASKER_comparativ_FC10_$tab1" . "_Regions_merged.tab".".gff", $threads);
	#kmasker::functions::add_annotation($fasta, "temp/KMASKER_comparativ_FC10_$tab2" . "_Regions_merged.tab", $BLASTDB, "temp/KMASKER_comparativ_FC10_$tab2" . "_Regions_merged.tab".".gff", $threads);
    system("cp" . " $temp_path/KMASKER_comparativ_FC10_$tab1" . "_Regions_merged".".gff" . " KMASKER_comparativ_FC10_$tab1" . ".gff");
    system("cp" . " $temp_path/KMASKER_comparativ_FC10_$tab2" . "_Regions_merged".".gff" . " KMASKER_comparativ_FC10_$tab2" . ".gff");
    #With Annotation
    #system("cp" . " temp/KMASKER_comparativ_FC10_$tab1" . "_Regions_merged_annotation".".gff" . " KMASKER_comparativ_FC10_$tab1" . "_annotation.gff");
    #system("cp" . " temp/KMASKER_comparativ_FC10_$tab2" . "_Regions_merged_annotation".".gff" . " KMASKER_comparativ_FC10_$tab2" . "_annotation.gff");

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

sub show_version_PM_run{
	print "\n\nVERSION of module run: ".$version_PM_run."\n\n";
}

	
	

1;
