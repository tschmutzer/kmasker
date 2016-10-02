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
		#my $href = &get_kindex_info($href_info, $href_repo);
		my %HASH_info_this = %{$href_info};
		
		# GET info
		my $rept = $HASH_info_this{"rept"};
		
		my @ARRAY_kindex_info 	= split("\t", $HASH_repo_this{$kindex});
		my $seq_depth			= $ARRAY_kindex_info[4];
		my $k					= $ARRAY_kindex_info[5];
		my $md5sum				= $ARRAY_kindex_info[8]; 
		my $absolut_path		= $ARRAY_kindex_info[10];		
		
		#create symbolic link to kindex
		system("ln -s ".$absolut_path."KINDEX_".$kindex."_k".$k.".jf"); 
		#start
		#system("cmasker -f ".$fasta." -j ".$kindex." -n ".$seq_depth." -r ".$rept." -o Kmasker_KDX_".$kindex."_RT".$rept);
		system("cmasker -f ".$fasta." -j KINDEX_".$kindex."_k".$k.".jf -n ".$seq_depth." -r ".$rept." -o");
		
		#clean
		system("rm KINDEX_".$kindex."_k".$k.".jf");
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