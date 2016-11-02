#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::occ;
my $occ;
my $repT;
my $fasta;

 GetOptions ("occ=s"   => \$occ,  
             "repT=i"  => \$repT,
             "fasta=s" => \$fasta) 
 or die("Error in command line arguments\n");

apply_occ($fasta, $occ, $repT)