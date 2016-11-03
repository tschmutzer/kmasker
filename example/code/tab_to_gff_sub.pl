#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::filehandler;
my $tab;
my $subtab;
my $feature;
my $min = 10;

 GetOptions ("tab=s"   => \$tab,
 			 "subtab=s" => \$subtab,  
             "feature=s" => \$feature) 
 or die("Error in command line arguments\n");

kmasker::filehandler::tab_to_gff($feature, $min, $tab, $subtab);