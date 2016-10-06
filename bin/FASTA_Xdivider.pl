#!/usr/bin/perl -w
use strict;
use IO::File;
use File::Basename;
use Getopt::Long;
use kmasker::filehandler;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author:       Thomas Schmutzer
# date:         2011_06_17
# last update:	2016_10_05 by Chris Ulpinnis
# institute:    @IPK Gatersleben
# version:		0.0.102
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# variables
my $fasta 		= "";
my $sizelimit 	= 20;
my $sl;
my $help;


my $result = GetOptions ("fasta=s"   	=> \$fasta,  		# provide the fasta file
						 "sl=s"			=> \$sl,
						 "help" 		=> \$help);


if(defined $help){
	print "\n\n USAGE of the programm FASTA_Xdivider.pl:";
	print "\n\t FASTA_Xdivider.pl --fasta sequences.fasta ";
	print "\n\n OPTIONS:\n"; 
	print "\n\t --fasta\t provide the fasta file";
	print "\n\n";
	exit();
}





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++ main ++++++++++++++++++++++++++++++++++++++


print "\n\t.. dividing the FASTA!!";
print "";


$sizelimit = $sl if(defined $sl);
&process_fasta;
print "\n\t.. finished\n";

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 



################################
# subroutin 
sub process_fasta(){
	
	
	# Initiating Handler	
	open( my $inFASTA, "<", "$fasta");
	(my $name,my $path,my $suffix) = fileparse($fasta, qr/\.[^.]*/);
	open( my $newFAST, ">", $path . "/Xsplit_" . $name . $suffix);
	my %seqdata; 	
	while(read_sequence($inFASTA, \%seqdata)) {
		
		my $id 			= $seqdata{header};
		my $seq			= $seqdata{seq};
		my @seqarray	= split("X", $seq);
		my $split 		= 1;
		
		foreach my $element (@seqarray){
		
			if(length($element) >= $sizelimit){
				my $newID = $id."_".$split++;
				my $desc  = length($element);
                print $newFAST ">".$id . " " . $desc ."\n";
                print $newFAST $element . "\n";
			}
		}	
        				
	}	
		   
   	close($inFASTA);
	close($newFAST);  	
	
}


exit(1);

