#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::filehandler;
my $tab;
my $feature;

 GetOptions ("tab=s"   => \$tab,  
             "feature=s" => \$feature) 
 or die("Error in command line arguments\n");

kmasker::filehandler::tab_to_gff($feature, $tab);