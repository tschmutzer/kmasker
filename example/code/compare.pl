#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::filehandler;
use kmasker::occ;

my $occ1;
my $occ2;
my $threshold;
my $rt = 40;


 GetOptions ("1=s"   => \$occ1,
 			 "2=s" => \$occ2,
 			 "t=i" => \$threshold,
 			 "r=i" => \$rt)
 or die("Error in command line arguments\n");

kmasker::occ::multi_occ($rt, $threshold, $occ1, $occ2, "test_");

