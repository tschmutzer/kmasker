#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::filehandler;
my $fasta;
my $prefix;

 GetOptions ("fasta=s"   => \$fasta,  
             "prefix=s" => \$prefix) 
 or die("Error in command line arguments\n");

kmasker::filehandler::fasta_to_tab($fasta, $prefix);