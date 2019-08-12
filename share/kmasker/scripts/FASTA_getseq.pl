#!/usr/bin/perl -w
use strict;
use IO::File;
use File::Basename;
use Getopt::Long;
use POSIX;


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author:       Thomas Schmutzer
# date:         2017_08_17
# last update:	2018_01_12
# institute:    @MLU
# version:		1.0.14
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
my $version = "v1.0.14 rc180112";

# variables
my $fasta;
my $occ;
my $help 	= "";
my $list;
my %HASH_LIST = ();

my $result = GetOptions ("fasta=s"   	=> \$fasta,  		# provide the fasta file
						 "occ=s"		=> \$occ,			# provided occ freq file
						 "list=s"		=> \$list,
						 "help" 		=> \$help);


if($help ne ""){
	print "\n\n _________________________________________________________";
	print "\n ---------------------------------------------------------";
	print "\n FASTA_getseq (implemented by Thomas Schmutzer):\n";
	print "\n version = ".$version."\n";
	print " ---------------------------------------------------------";
	print "\n\n USAGE of the programm FASTA_getseq.pl:";
	print "\n\t FASTA_getseq.pl --fasta sequences.fasta --rlist IDfile.txt\n";
	print "\n\n OPTIONS:\n"; 
	print "\n\t --fasta-\t provide the fasta file";
	print "\n\t --occ\t-\t provided occ file";
	print "\n\t --list-\t provide file of sequence identifier (one ID per line)";
	exit();
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++ main ++++++++++++++++++++++++++++++++++++++


print "\n\t .. start extraction of sequence!!\n";

&read_IDLIST();
&process();

print "\n\t .. selection process finished\n\n";

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
################################
# subroutin 
sub process(){
	
	# Initiating Handler
	my $file = ""; 
	$file = $occ if(defined $occ);
	$file = $fasta if(defined $fasta);
	my $FILE_IN		= new IO::File($file, "r") or die "unable to read ".$file." \n\n";	
	my $FILE_PRINT 	= new IO::File($file.".selection", "w") or die "unable to read ".$file.".selection \n\n";
		print "\n\t .. extracting list of sequences from ".$file."\n";
	
	my $this_id = "";
	my $status 	= 0;	# 1 = write; 0 = do not write 
	while(<$FILE_IN>){
		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		my $line = $_;
		my $original_line = $line;
		$line =~ s/\n$//;
	
		#processing
		if($line =~ /^>/){
			#id			
			my @array = split(" ", $line);
			$this_id = $array[0];
			$this_id =~ s/^>//;
			if(exists $HASH_LIST{$this_id}){
				$status = 1;
			}else{
				$status = 0;
			}
		}
		
		if($status == 1){
			print $FILE_PRINT $original_line;
		}
	}		
				   
   	$FILE_IN->close();
   	$FILE_PRINT->close();
   	
}



#################################
# subroutin should read in the listfile and collect the readID with 
# corresponding start\tend positions in the RIDS Hash
#
sub read_IDLIST(){
	my $LIST = new IO::File($list, "r") or die "\n unable to read $list file $!";
	while(<$LIST>){
		my $line = $_;
		next if($line =~ /^$/);
		$line =~ s/\n|\t|\s//g;
		$HASH_LIST{$line} = 1;		
	}
}


