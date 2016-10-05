package kmasker::kmasker_postprocessing;
use Exporter qw(import);
use File::Basename;
use strict;
use warnings;

our @ISA = qw(Exporter);
our @EXPORT = qw(
plot_histogram
repeat_annotation
gff_construction
);
our @EXPORT_OK = qw(plot_histogram repeat_annotation gff_construction);


## subroutine
#
sub plot_histogram{
    $occ=$_[0];
    $clist=$_[1];
    system("occVisualizer.R -i ".$occ." -l ".$clist);
    #The script will skip large contigs to avoid long running times
    #You can force it to do it anyway with -f
    #If you do not want to provide a contig list (instead running on all entries
    #in the occ file) use the following
    #system("occVisualizer.R -i ".$occ.");
    #The script will generate one plot per contig in the current working directory
}


## subroutine
#
sub gff_construction{
#implement
	
}

## subroutine
#
sub repeat_annotation{
#implement	

}




1;