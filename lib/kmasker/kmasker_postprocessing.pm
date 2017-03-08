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

## VERSION
my $version_PM_postprocessing 	= "0.0.1 rc170308";

## subroutine
#
sub plot_histogram{
    my $occ		=	$_[0];
    my $clist	=	$_[1];
    system("occVisualizer.R -i ".$occ." -l ".$clist);
    #The script will skip large contigs to avoid long running times
    #You can force it to do it anyway with -f
    #If you do not want to provide a contig list (instead running on all entries
    #in the occ file) use the following
    #system("occVisualizer.R -i ".$occ.");
    #The script will generate one plot per contig in the current working directory
    
    #
    my @ARRAY_info = split("_", $occ);
    my $kindex = $ARRAY_info[1];
    if(!( -d "./Kmasker_plots_".$kindex)){
    	system("mkdir Kmasker_plots_".$kindex);
    }
   
   	my $LIST = new IO::File($clist, "r") or die "\n unable to read $clist $!";	
   	while(<$LIST>){
   		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		my $line = $_;
		$line =~ s/\n//;
		system("mv ".$line."*.png Kmasker_plots_".$kindex);
   	}    
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