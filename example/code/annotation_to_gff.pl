#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::filehandler;
my $gff;
my $blast;

 GetOptions ("gff=s"   => \$gff,  
             "blast=s" => \$blast) 
 or die("Error in command line arguments\n");

kmasker::filehandler::add_annotation_to_gff($gff, $blast);