#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::occ;
my $tab;
my $percent;
my $min = 5;
 GetOptions ("tab=s"   => \$tab,  
             "percent=i" => \$percent) 
 or die("Error in command line arguments\n");

kmasker::filehandler::merge_tab_seeds($tab, $percent, $min);