package kmasker::functions;
use kmasker::occ;
use kmasker::filehandler;
use File::Basename;
use Exporter qw(import);
use strict;
use warnings;
use Data::Uniqid qw ( suniqid luniqid );
our @ISA = qw(Exporter);
our @EXPORT = qw(fasta_to_uppercase Xtract add_annotation getLoggingTime getsID getlID);
our @EXPORT_OK = qw(fasta_to_uppercase Xtract add_annotation getLoggingTime getsID getlID);


## VERSION
my $version_PM_functions 	= "0.0.2 rc181025";

#sub add_repository {
   # function adds local repository to global shared repository 
#}


#sub show_repository {
   # return list of available k-mer indices
   
   
#}

#sub read_write_repository {
	#read config 
	
#	my $usr = `echo \$USER`;
#	print "\n USER is: ".$usr;
	
#}

sub getsID {
	return(suniqid);
}
sub getlID {
	return(luniqid);
}

sub getLoggingTime {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d_%02d%02d%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}

sub fasta_to_uppercase {
	my $fasta 		= $_[0];
	open( my $inFASTA, "<", "$fasta");
	(my $name,my $path,my $suffix) = fileparse($fasta, qr/\.[^.]*/);
	open( my $newFASTA, ">", $path . "/UC_" . $name . $suffix);
	my %seqdata; 	
	while(read_sequence($inFASTA, \%seqdata)) {
		
		my $id 			= $seqdata{header};
		my $seq			= $seqdata{seq};
		print $newFASTA ">$id\n" . uc($seq) ."\n";

   				
	}	
		   
   	close($inFASTA);
	close($newFASTA);
	return( $path . "/UC_" . $name . $suffix);

}

sub add_annotation {
        my $FASTA = $_[0]; #FASTA to extract sequences from
        my $BLAST_db = $_[1]; #BLAST-reference
        my $GFF = $_[2]; #GFF-file to be annotated
        my $feature = $_[3];
        my %HASH_info = %{$_[4]};
        my $threads = $HASH_info{"threads"};
        my $temp_dir = $HASH_info{"temp_path"};
        kmasker::filehandler::extract_feature_gff($FASTA, $GFF, $feature, $temp_dir);
        #extract_sequence_region($FASTA, $TAB);
        #system("mv selected_* temp/");
        # Using standard word size of megablast [28]
        if(exists $HASH_info{"user setting blast"}) {
       	      my $parameterstring = $HASH_info{"user setting blast"};
       	      system("blastn -db \"" . $BLAST_db . "\" -query " . "${temp_dir}/selected_" . $FASTA . " -num_threads ".$threads." -outfmt 6 " . $parameterstring . " -ungapped -max_hsps 1 -max_target_seqs 1" . "  -out ${temp_dir}/kmasker_blast.txt");
        }
        else{
        	system("blastn -db \"" . $BLAST_db . "\" -query " . "${temp_dir}/selected_" . $FASTA . " -perc_identity 80 -evalue 0.1 -num_threads ".$threads." -outfmt 6 -ungapped -max_hsps 1 -max_target_seqs 1" . " -out ${temp_dir}/kmasker_blast.txt");
        }
        #FIXME: Add exchangeable configuration to blast
        #--->almost done, just add the necessary parameters to %HASH_INFO
        kmasker::filehandler::add_annotation_to_gff($GFF, "${temp_dir}/kmasker_blast.txt");
}

sub Xtract{
	
	my $fasta 		= $_[0];
	my $sizelimit 	= $_[1];
	my $suffix		= $_[2];
	my $N_ration 	= $_[3];
	
	if(!defined $N_ration){
		$N_ration = 0.05;
	}
	
	# Initiating Handler	
	open( my $inFASTA, "<", "$fasta");
	(my $name,my $path,my $fsuffix) = fileparse($fasta, qr/\.[^.]*/);
	open( my $newFAST, ">", $path . "/KMASKER_filtered_regions_" . $suffix . ".fasta");
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
				my $str = $element;
				my $nr_N = $str =~ tr/N//;
				if(($nr_N/$desc)<$N_ration){
                	print $newFAST ">".$newID . " " . $desc ."\n";
                	print $newFAST $element . "\n";
				}
			}
		}       				
	}	
		   
   	close($inFASTA);
	close($newFAST);  	
	
}
