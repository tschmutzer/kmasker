package kmasker::kmasker_explore;
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
our @EXPORT_OK = qw(plot_histogram custom_annotation gff_construction);

## VERSION
my $version_PM_explore 	= "0.0.1 rc180419";

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

## subroutine
#
sub custom_annotation{
	my $fasta		=	$_[0];
    my $gff			=	$_[1];
    my $href_DB 	=   $_[2]; 
    my $href_info	=	$_[3];
    
    my %HASH_info 	= %{$href_info};
    my %HASH_DB		= %{$href_DB};
    my $db_fasta;
    my $db;
    
    if(exists $HASH_DB{"db"}){
    	# check if file is a blastable DB
        if((! -e $db . ".nhr"  ) || (! -e $db .".nin") || (! -e $db . ".nsq")) {
            print("The BLASTdb was not found or is not complete!\n");
            print("Please run Kmasker with dbfasta instead of db again to rebuild the BLASTdb!\n");   
        } 
        $db = $HASH_DB{"db"};    	  	
    }
    elsif(exists $HASH_DB{"db_fasta"}){
        $db_fasta = $HASH_DB{"db_fasta"};
        $db=$db_fasta;
         if((! -e $db . ".nhr"  ) || (! -e $db .".nin") || (! -e $db . ".nsq")) {
            print("BLASTdb is missing. It will be built now!\n");
            system("makeblastdb -in \"".$db."\" -dbtype nucl ");
            print("BLASTdb was built. You can use the path to your fasta just with -db in the future.\n");   
        }
    }
    else{
    	print "\n\n Parameter dbfasta or db have to be provided for custom annotation!";
    	print   "\n Kmasker was stopped.\n\n";
    	exit();
    }
    
    
    # TASK for Chris to move annotation from run to this section!
  	my $threads		= $HASH_info{"threads"};	
    my $temp_path 	= $HASH_info{"temp_path"};   
	kmasker::functions::add_annotation($fasta, ".tab", $db, "$temp_path/.gff", $threads);
    #ToDo Reimplemantaion of blast function 
    #Use of extract_feature_gff (feature must be a feature given by the user)
    
	
}





1;