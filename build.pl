#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::occ;
my $fasta;
my $index;
my $mer;

 GetOptions ("k=i" => \$mer,    # numeric
             "fasta=s"   => \$fasta,      # string
             "index=s"  => \$index)   # flag
 or die("Error in command line arguments\n");

kmasker::occ::make_occ($fasta, $index, $mer)