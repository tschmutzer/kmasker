package kmasker::functions;
use Exporter qw(import);
use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_write_repository);
our @EXPORT_OK = qw(read_write_repository);


## VERSION
my $version_PM_functions 	= "0.0.1 rc170308";

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
