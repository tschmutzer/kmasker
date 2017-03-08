package kmasker::kmasker_run;
use Exporter qw(import);
use File::Basename;
use strict;
use warnings;
use kmasker::filehandler;
use kmasker::occ;

#adapt
our @ISA = qw(Exporter);
our @EXPORT = qw(
	run_kmasker_SK
	run_kmasker_MK
	show_version_PM_run
);
our @EXPORT_OK = qw(run_kmasker_SK run_kmasker_MK show_version_PM_run);


## VERSION
my $version_PM_run 	= "0.0.24 rc170308";

## subroutine
#
sub run_kmasker_SK{
	# SINGLE KINDEX (SK)
	my $fasta 		= $_[0];
	my $kindex		= $_[1];
	my $href_info	= $_[2];
	my $href_repo 	= $_[3];
	my %HASH_repo_this = %{$href_repo};
   	
	if(exists $HASH_repo_this{$kindex}){
		
		# GET info from repository
		my %HASH_info_this 		= %{$href_info};
		
		# GET info
		my $rept 				= $HASH_info_this{"rept"};
		my $length_threshold 	= $HASH_info_this{"min_length"};
		
		my @ARRAY_kindex_info 	= split("\t", $HASH_repo_this{$kindex});
		my $seq_depth			= $ARRAY_kindex_info[4];
		my $k					= $ARRAY_kindex_info[5];
		my $md5sum				= $ARRAY_kindex_info[8]; 
		my $absolut_path		= $ARRAY_kindex_info[10];		
		
		#create symbolic link to kindex from private or global
		my $full_kindex_name = "KINDEX_".$kindex."_".$md5sum."_k".$k.".jf";		
		if(-e $absolut_path."/".$full_kindex_name){
			system("ln -s ".$absolut_path."/".$full_kindex_name);
		}else{
			print "\n WARNING: KINDEX not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
			exit();
		}
		
		#start
		system("cmasker -f ".$fasta." -j ".$full_kindex_name." -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex);
       
       #clean
		system("rm ".$full_kindex_name);
        if(!(-e "KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta)) {
            print "\n .. KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta." was not generated!\n";
            print "\n .. please provide a bug report!\n\n";
            exit();
        }
        system("FASTA_Xdivider.pl --fasta KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta." --sl ".$length_threshold);
        system("mv Xsplit_KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta." KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_MIN".$length_threshold."_".$fasta);
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
	my $href_info		= $_[2];
	my $href_repo 		= $_[3];
	my %HASH_repo_this 	= %{$href_repo};
	my @ARRAY_kindex 	= @{$aref_mkindex};
	
	# GET info from repository
	my %HASH_info_this 		= %{$href_info};
		
	# GET info
	my $rept 				= $HASH_info_this{"rept"};
	my $length_threshold 	= $HASH_info_this{"min_length"};
	my @ARRAY_full_kindex_names = ();
	my @ARRAY_seq_depth			= ();
	
	#ADD info
	my @ARRAY_occ = ();
	
	#USER info
	my $fold_change = 10;
	
	
	foreach(@ARRAY_kindex){
   		my $kindex = $_;
		if(exists $HASH_repo_this{$kindex}){		
			
			my @ARRAY_kindex_info 	= split("\t", $HASH_repo_this{$kindex});
			my $seq_depth			= $ARRAY_kindex_info[4];
			my $k					= $ARRAY_kindex_info[5];
			my $md5sum				= $ARRAY_kindex_info[8]; 
			my $absolut_path		= $ARRAY_kindex_info[10];		
			
			#create symbolic link to kindex from private or global
			my $full_kindex_name = "KINDEX_".$kindex."_".$md5sum."_k".$k.".jf";
			push(@ARRAY_full_kindex_names, $full_kindex_name);
			push(@ARRAY_seq_depth, $seq_depth);
			if(-e $absolut_path."/".$full_kindex_name){
				system("ln -s ".$absolut_path."/".$full_kindex_name);
				
				#PRODUCE OCC
				system("cmasker -f ".$fasta." -j ".$full_kindex_name." -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex." -s");
				
				my @NAME = split(/\./, $fasta);
				pop @NAME;
				my $name = join(".", @NAME);
				print "\n NAME = ".$name;
				push(@ARRAY_occ, "KMASKER_".$kindex."_N".$seq_depth."_".$name.".occ");
				
				#clean
    			system("rm ".$full_kindex_name);
				
			}else{
				print "\n WARNING: KINDEX not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
				exit();
			}
		}else{
			#KINDEX is missing in repository
			print "\n .. Kmasker was stopped!\n";
			print "\n .. The kindex (".$kindex.") you requested is not available in repository!\n\n";
		}	
	}
		
	#Convertion to TAB file
	my $occ1 = $ARRAY_occ[0];
	my $occ2 = $ARRAY_occ[1];
	kmasker::occ::multi_occ($rept, $fold_change, $occ1, $occ2, "KMASKER_comparativ_FC".$fold_change."_");
	(my $name1,my $path1,my $suffix1) = fileparse($occ1, qr/\.[^.]*/);
	(my $name2,my $path2,my $suffix2) = fileparse($occ2, qr/\.[^.]*/);
	
	foreach(@ARRAY_kindex){
   		my $kindex = $_;
		
		#EXTEND regions
		my $percent = 10; 	# 10% of length 	FIXME: That parameter has to come from user
		my $min_seed= 5;	# 5 bp 				FIXME: That parameter has to come from user
		kmasker::filehandler::merge_tab_seeds($name1.".tab", $percent, $min_seed);
		
		#PRODUCE GFF
		my $min_gff= 10;	# 10 bp minimal length to be reported in GFF
		my $feature = "MDR";
		kmasker::filehandler::tab_to_gff($feature, $min_gff, $name1.".tab_growed");
		
	}
	
	     
   #clean

#    if(!(-e "KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta)) {
#    	print "\n .. KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta." was not generated!\n";
#        print "\n .. please provide a bug report!\n\n";
#        exit();
#    }
#       system("FASTA_Xdivider.pl --fasta KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta." --sl ".$length_threshold);
#        system("mv Xsplit_KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_".$fasta." KMASKER_".$kindex."_RT".$rept."_N".$seq_depth."_MIN".$length_threshold."_".$fasta);
#	}
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