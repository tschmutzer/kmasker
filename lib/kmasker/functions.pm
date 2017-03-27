package kmasker::functions;
use kmasker::occ;
use kmasker::filehandler;
use File::Basename;
use Exporter qw(import);
use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_write_repository Xtract add_annotation);
our @EXPORT_OK = qw(read_write_repository Xtract add_annotation);


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

sub add_annotation {
        my $FASTA = $_[0]; #FASTA to extract sequences from
        my $TAB = $_[1]; #TAB-file with regions to extract
        my $BLAST_db = $_[2]; #BLAST-reference
        my $GFF = $_[3]; #GFF-file to be annotated
        #FIXME: Give configuration hash as 5th parameter.

        extract_sequence_region($FASTA, $TAB);
        system("mv selected_* temp/");
        system("blastn -db " . $BLAST_db . " -query " . "temp/selected_" . $FASTA . " -perc_identity 80 -word_size 50 -evalue 0.1 -num_threads 30 -outfmt 6 -ungapped -max_hsps 1 -max_target_seqs 1" . " -out temp/kmasker_blast.txt");
        #FIXME: Add exchangeable configuration to blast
        add_annotation_to_gff($GFF, "temp/kmasker_blast.txt");
}

sub Xtract{
	
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
