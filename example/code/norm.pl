#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::occ;
my $occ;
my $depth;

 GetOptions ("occ=s"   => \$occ,  
             "depth=i"  => \$depth) 
 or die("Error in command line arguments\n");

normalize_occ($occ, $depth)