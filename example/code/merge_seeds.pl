#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::occ;
my $occ1;
my $percent;

 GetOptions ("occ1=s"   => \$occ1,  
             "percent=i" => \$percent) 
 or die("Error in command line arguments\n");

kmasker::occ::merge_seeds($occ1, $percent);