#!/usr/bin/perl -w
use strict;
use IO::File;
use File::Basename;
use Bio::SeqIO;
use Getopt::Long;



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: 	Thomas Schmutzer
# Dev. date	2010_02_24
# UPDATE:	2016_05_03
# institute: @IPK Gatersleben
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


my ($fasta) = ("");
my $help = "";
my $list;
my $out	 = "";



my $result = GetOptions ("fasta=s"	=>\$fasta,
                         "list=s"	=>\$list,
                         "out"		=>\$out,		# set out if a file should be created containing sequences which are not empty (seq.length  != 0)
                         "help" 	=>\$help
                         );    

if($help ne ""){
   &info4help
}

my $fasta_in 	= Bio::SeqIO->new(-file => $fasta , '-format' => 'fasta');
my $fasta_out 	= Bio::SeqIO->new(-file => ">non0_".$fasta , '-format' => 'fasta') or die "unable to write non0_".$fasta;
my $LENGTH 		= new IO::File($fasta.".length", "w") or die "unable to writre!\n";


#####
# reading selection list
my %SEL;
if(defined $list){
	my $list_in 	= new IO::File($list, "r") or die "unable to read $list !\n";
	while(<$list_in>){
		my $line = $_;
		next if(($line =~ /^$/)||($line =~ /^#/));
		$line =~ s/\n//;
		$SEL{$line} = 1;
	}
}
   	
######
# compute length by using BioPerl   	
while (my $read = $fasta_in->next_seq()){
   		
   	my $ID 		= $read->id();
   	if(defined $list){
   		if(exists $SEL{$ID}){
   			my $length 	= $read->length();  
	   		if($length != 0){
	   			$fasta_out->write_seq($read) if($out ne "");
	   		}	   	
	   		print $LENGTH $ID."\t".$length."\n";	
   		}   		
   	}else{
	   	my $length 	= $read->length();  
	   	if($length != 0){
	   		$fasta_out->write_seq($read) if($out ne "");
	   	}	   	
	    print $LENGTH $ID."\t".$length."\n";
   	}
}
   
$fasta_in->close();
system("rm non0_".$fasta); 
 
 
#=============================================================================
sub info4help() {

my $prog = basename($0); 
print STDERR <<EOF ;

 ======= $prog ========
#
# Description:

perl $prog --help
Note:
        options:  
        
        --fasta	FASTA file with reads,
        --list	list of IDs for selective mode
        (only reads in the list are computed)


EOF
exit(1) ;

}

1;

#END




