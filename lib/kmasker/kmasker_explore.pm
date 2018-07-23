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
our @EXPORT_OK = qw(plot_histogram_raw plot_histogram_mean custom_annotation gff_construction report_statistics);

## VERSION
my $version_PM_explore 	= "0.0.1 rc180720";

## subroutine
#
sub plot_histogram_mean{
    my $occ		=	$_[0];
    my $list	=	$_[1];
    
    #VAR
    my $outlist = "kmasker_seq.ids";
    
    if(defined $list){
    #SUBSET
	    my $systemcall_wc = `wc -l $list`;  
	    $systemcall_wc =~ s/ /\t/;
	    $systemcall_wc =~ s/ //g;
	  	my @ARRAY_sys = split("\t", $systemcall_wc);
		my $number = $ARRAY_sys[0];
		print "\n ".$number." sequences selected for plotting.";
		if($number >= 10000){
			print "\n WARNING: Generating more than 10000 plots will take long and is not reccomended!\n";
		}   
	
		#subset
		system("FASTA_getseq.pl -o ".$occ." --list ".$list);
		my $occ_subset = $occ.".selection";
    
	    #vis
	    system("occVisualizer.R -i ".$occ_subset." -l ".$list);
	    #The script will generate one plot per contig in the current working directory
	    #The script will skip large contigs to avoid long running times
	    #You can force it to do it anyway with -f
	    
	    #clean
   		system("rm ".$occ_subset);
   		
    }else{
    #FULL DATASET
    	$list = $outlist;
    	system("grep \">\" ".$occ."| sed \'s/^>//\' | awk -F\" \" '{print \$1}' >".$outlist);
    	system("occVisualizer.R -i ".$occ);
    	#The script will generate one plot per contig in the current working directory
	    #If you do not want to provide a contig list (instead running on all entries
	    #in the occ file) use the following
	    
    }
    
    #STORE results
    my @ARRAY_info = split("_", $occ);
    my $kindex = $ARRAY_info[1];
    if(!( -d "./Kmasker_plots_".$kindex)){
    	system("mkdir Kmasker_plots_".$kindex);
    }   
   	my $LIST = new IO::File($list, "r") or die "\n unable to read $list $!";	
   	while(<$LIST>){
   		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		my $line = $_;
		$line =~ s/\n//;
		if(-e $line.".png"){
			system("mv ".$line.".png Kmasker_plots_".$kindex."/".$line."_hist.png");
		}
   	}
   	
   	#clean    
   	if(-e "kmasker_seq.ids"){
 		system("rm kmasker_seq.ids");
   	}
}


## subroutine
#
sub plot_histogram_raw{
    my $occ		=	$_[0];
    my $list	=	$_[1];    
    
}


## subroutine
#
sub report_statistics{
#implement
	
}


## subroutine
#
sub custom_annotation{
	my $fasta		=	$_[0];
    my $gff			=	$_[1];
    my $feature     =   $_[2];
    my $href_DB 	=   $_[3]; 
    my $href_info	=	$_[4];
    
    my %HASH_info 	= %{$href_info};
    my %HASH_DB		= %{$href_DB};
    my $db_fasta;
    my $db;
    mkdir($HASH_info{"temp_path"}, 0755);
    
    if(exists $HASH_DB{"db"}){
    	# check if file is a blastable DB
        $db=$HASH_DB{"db"};
        if((! -e $db . ".nhr"  ) || (! -e $db .".nin") || (! -e $db . ".nsq")) {
            print("The BLASTdb was not found or is not complete!\n");
            print("Please run Kmasker with dbfasta instead of db again to rebuild the BLASTdb!\n");   
        } 
    }
    elsif(exists $HASH_DB{"db_fasta"}){
        $db_fasta = $HASH_DB{"db_fasta"};
        my($db_prefix, $dirs, $suffix) = fileparse($db_fasta, (".fa", ".fasta"));
        $db=$dirs."/".$db_prefix.$suffix;
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
    
    kmasker::functions::add_annotation($fasta, $db, $gff, $feature ,$href_info);
    (my $name,my $path,my $suffix) = fileparse($gff, qr/\.[^.]*/);
    if(-x "$path/${name}_with_annotation${suffix}") {
      print("\nAnnotation of GFF was done!\n");
      unlink($gff);
      system("mv" . " " . "$path/${name}_with_annotation${suffix}" . " " . $gff);
    }
    else{
        print("An annotated GFF was not created. Something went wrong!");
    }
	
}





1;
