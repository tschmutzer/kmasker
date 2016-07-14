#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use kmasker::filehandler;
my $occ;

 GetOptions ("occ=s"   => \$occ)      # string
 or die("Error in command line arguments\n");

occ_length($occ);