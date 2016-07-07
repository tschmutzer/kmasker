#!/usr/bin/perl -w
use strict;
use IO::File;
use File::Basename;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author:       Thomas Schmutzer
# date:         2011_06_17
# last update:	2011_06_17
# institute:    @IPK Gatersleben
# version:		0.0.9
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
	my $inFASTA		= new Bio::SeqIO(-file =>$fasta, -format =>'fasta') or die "could not find file $fasta \n";
	my $newFAST		= new Bio::SeqIO(-file =>">Xsplit_".$fasta, -format =>'fasta') or die "could not find file ".$fasta.".qual\n";		
		
	while(my $read = $inFASTA->next_seq()) {
		
		my $id 			= $read->id();
		my $seq			= $read->seq();
		my @seqarray	= split("X", $seq);
		my $split 		= 1;
		
		foreach my $element (@seqarray){
		
			if(length($element) >= $sizelimit){
				my $newID = $id."_".$split++;
				my $desc  = length($element);
				my $seqobj = Bio::PrimarySeq->new(
						 '-seq' 	=> $element,
		                 '-id'  	=> $newID,
		                 '-desc' => $desc
		                );		                
		        $newFAST->write_seq($seqobj);
			}
		}	
        				
	}	
		   
   	$inFASTA->close();
   	$newFAST->close();   	
   	
	
}


exit(1);

