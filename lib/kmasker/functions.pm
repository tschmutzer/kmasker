package kmasker::functions;
use Exporter qw(import);
use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_write_repository Xtract);
our @EXPORT_OK = qw(read_write_repository Xtract);


## VERSION
my $version_PM_functions 	= "0.0.1 rc170310";

sub add_repository {
   # function adds local repository to global shared repository 
}


sub show_repository {
   # return list of available k-mer indices
   
   
}

sub read_write_repository {
	#read config 
	
	my $usr = `echo \$USER`;
	print "\n USER is: ".$usr;
	
}


################################
# subroutin 
sub Xtract(){
	
	my $fasta 		= $_[0];
	my $sizelimit 	= 20;	#FIXME HASH_info	
	
	# Initiating Handler	
	open( my $inFASTA, "<", "$fasta");
	(my $name,my $path,my $suffix) = fileparse($fasta, qr/\.[^.]*/);
	open( my $newFAST, ">", $path . "/Xsplit_" . $name . $suffix);
	my %seqdata; 	
	while(read_sequence($inFASTA, \%seqdata)) {
		
		my $id 			= $seqdata{header};
		my @ARRAY_id	= split(" ", $id);
		$id				= $ARRAY_id[0];	
		my $seq			= $seqdata{seq};
		my @seqarray	= split("X", $seq);
		my $split 		= 1;
		
		foreach my $element (@seqarray){
		
			if(length($element) >= $sizelimit){
				my $newID = $id."_".$split++;
				my $desc  = length($element);
                print $newFAST ">".$newID . " " . $desc ."\n";
                print $newFAST $element . "\n";
			}
		}       				
	}	
		   
   	close($inFASTA);
	close($newFAST);  	
	
}
