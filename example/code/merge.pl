#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::occ;
my $occ1;
my $occ2;
my $occout;

 GetOptions ("occ1=s"   => \$occ1,  
             "occ2=s"  => \$occ2,
             "occout=s" => \$occout) 
 or die("Error in command line arguments\n");

merge_occ($occ1, $occ2, $occout);