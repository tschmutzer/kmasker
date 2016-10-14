package kmasker::kmasker_run;
use Exporter qw(import);
use File::Basename;
use strict;
use warnings;

#adapt
our @ISA = qw(Exporter);
our @EXPORT = qw(
	run_kmasker
);
our @EXPORT_OK = qw(run_kmasker);


## subroutine
#
sub run_kmasker{
	my $fasta 		= $_[0];
	my $kindex		= $_[1];
	my $href_info	= $_[2];
	my $href_repo 	= $_[3];
	my %HASH_repo_this = %{$href_repo};
   	
	if(exists $HASH_repo_this{$kindex}){
		
		# GET info from repository
		my %HASH_info_this = %{$href_info};
		
		# GET info
		my $rept 				= $HASH_info_this{"rept"};
		my $length_threshold 	= $HASH_info_this{"min_length"};
		
		my @ARRAY_kindex_info 	= split("\t", $HASH_repo_this{$kindex});
		my $seq_depth			= $ARRAY_kindex_info[4];
		my $k					= $ARRAY_kindex_info[5];
		my $md5sum				= $ARRAY_kindex_info[8]; 
		my $absolut_path		= $ARRAY_kindex_info[10];		
		
		#create symbolic link to kindex from private
		#FIXME 
		if(-e $HASH_info_this{"PATH_kindex_private"}."KINDEX_".$kindex."/KINDEX_".$kindex."_".$md5sum."_k".$k.".jf"){
			system("ln -s ".$HASH_info_this{"PATH_kindex_private"}."KINDEX_".$kindex."/KINDEX_".$kindex."_".$md5sum."_k".$k.".jf");
		}elsif(-e $HASH_info_this{"PATH_kindex_global"}."KINDEX_".$kindex."/KINDEX_".$kindex."_".$md5sum."_k".$k.".jf"){
			system("ln -s ".$HASH_info_this{"PATH_kindex_global"}."KINDEX_".$kindex."/KINDEX_".$kindex."_".$md5sum."_k".$k.".jf");
		}else{
			print "\n WARNING: KINDEX not found in path. Please check path variables! \n\t Kmasker has been stopped\n\n";
		}
		#start
		system("cmasker -f ".$fasta." -j KINDEX_".$kindex."_".$md5sum."_k".$k.".jf -n ".$seq_depth." -r ".$rept." -o" . " -p" .$kindex);
        #Note: It is also possible to use the absolute path to the KINDEX
		#clean
#		system("rm KINDEX_".$kindex."_".$md5sum."_k".$k.".jf");
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
#
sub get_kindex_info{
	my $href_info_this = $_[0];
	my $href_repo_this = $_[1];
	my %HASH_info_tmp = %{$href_info_this};
	my %HASH_repo_tmp = %{$href_repo_this};
	### not required
	
}

	
	

1;