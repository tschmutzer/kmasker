#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::occ;
my $occ;
my $depth;
my $fasta;

 GetOptions ("occ=s"   => \$occ,  
             "depth=i"  => \$depth,
             "fasta=s" => \$fasta) 
 or die("Error in command line arguments\n");

apply_occ($fasta, $occ, $depth)