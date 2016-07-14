#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::filehandler;
my $fasta;

 GetOptions ("fasta=s"   => \$fasta)      # string
 or die("Error in command line arguments\n");

sequence_length($fasta)