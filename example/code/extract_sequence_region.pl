#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::filehandler;
my $fasta;
my $list;
my $offset = 0;

 GetOptions ("fasta=s"   => \$fasta,
 			 "list=s" => \$list,  
             "offset=i" => \$offset) 
 or die("Error in command line arguments\n");

kmasker::filehandler::extract_sequence_region($fasta, $list, $offset);