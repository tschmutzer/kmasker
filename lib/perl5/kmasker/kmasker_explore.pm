package kmasker::kmasker_explore;
use Exporter qw(import);
use File::Basename;
use File::Copy;
use strict;
use warnings;
use kmasker::filehandler;
use kmasker::occ;
use kmasker::functions;
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use POSIX; 

my $timestamp = getLoggingTime();
our $log = "log_explore_" . $kmasker::functions::PID . ".txt";
our $PID = $kmasker::functions::PID;
our @ISA = qw(Exporter);
our @EXPORT = qw(
    plot_histogram
    repeat_annotation
    gff_construction
    $log
    $PID
);
our @EXPORT_OK = qw(plot_histogram_raw plot_histogram_mean custom_annotation gff_construction report_statistics plot_maker plot_maker_direct plot_barplot);

## VERSION
my $version_PM_explore 	= "0.0.2 rc180815";

## DATE
my $loctime = localtime;
#$loctime = strftime('%Y%m%d_%H%M',localtime); ## outputs 1208171008
$loctime = $PID;

my $path        = dirname abs_path $0;      

## subroutine
#
sub plot_histogram_mean{
    my $occ		=	$_[0];
    my $list	=	$_[1];
    my $dynamic =   $_[2];
    my $force   =   $_[3];
    my $sws     =   $_[4];
    #VAR
    my $arguments = "";
    if(defined $dynamic) {
        $arguments = $arguments . " -d";
    }
    if(defined $force) {
        $arguments = $arguments . " -f";
    }
    if(defined $sws) {
        $arguments = $arguments . " -w $sws";
    }
    if(defined $log) {
        $arguments = $arguments . " -g";
    }
    my $outlist = "kmasker_seq_$PID.ids";
    
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
		system($path . "/" ."FASTA_getseq.pl -o ".$occ." --list ".$list. " >>$log 2>&1");
		my $occ_subset = $occ.".selection";
    
	    #vis
	    system($path . "/" ."occVisualizer.R -i ".$occ_subset." -l ".$list . "$arguments". " >>$log 2>&1");
	    #The script will generate one plot per contig in the current working directory
	    #The script will skip large contigs to avoid long running times
	    #You can force it to do it anyway with -f
	    
	    #clean
   		system("rm ".$occ_subset);
   		
    }else{
    #FULL DATASET
    	$list = $outlist;
    	system("grep \">\" ".$occ."| sed \'s/^>//\' | awk -F\" \" '{print \$1}' >".$outlist);
    	system($path . "/" ."occVisualizer.R -i ".$occ . "$arguments". " >>$log 2>&1");
    	#The script will generate one plot per contig in the current working directory
	    #If you do not want to provide a contig list (instead running on all entries
	    #in the occ file) use the following
	    
    }
    
    #STORE results
    if(!( -d "./kmasker_plots_$PID")){
    	mkdir("kmasker_plots_$PID");
    }   
   	my $LIST = new IO::File($list, "r") or die "\n unable to read $list $!";	
   	while(<$LIST>){
   		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		my $line = $_;
		$line =~ s/\n//;
		if(-e $line.".png"){
			system("mv ".$line.".png kmasker_plots_$PID". "/".$line."_histm.png");
		}
   	}
   	
   	#clean    
   	if(-e "kmasker_seq_$PID.ids"){
 		system("rm kmasker_seq_$PID.ids");
   	}
}


## subroutine
#
sub plot_histogram_raw{
    my $occ		=	$_[0];
    my $list	=	$_[1];    
    my $force   =   $_[2];
    #VAR
    my $arguments = "";
    if(defined $force) {
        $arguments = $arguments . " --force";
    }
    
    my $outlist = "kmasker_seq_$PID.ids";
    
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
        system($path . "/" ."FASTA_getseq.pl -o ".$occ." --list ".$list. " >>$log 2>&1");
        my $occ_subset = $occ.".selection";
    
        #vis
        my $call = $path . "/" ."occHistolizer.R -i ".$occ_subset." -l ".$list . "$arguments". " >>$log 2>&1";
        #print "CALL: " .$call;
        system($call);
        #The script will generate one plot per contig in the current working directory
        #The script will skip large contigs to avoid long running times
        #You can force it to do it anyway with -f
        
        #clean
        system("rm ".$occ_subset);
        
    }else{
    #FULL DATASET
        $list = $outlist;
        system("grep \">\" ".$occ."| sed \'s/^>//\' | awk -F\" \" '{print \$1}' >".$outlist);
        my $call = $path . "/" ."occHistolizer.R -i ".$occ . "$arguments". " >>$log 2>&1";
        #print "CALL: " .$call;
        system($call);
        #The script will generate one plot per contig in the current working directory
        #If you do not want to provide a contig list (instead running on all entries
        #in the occ file) use the following
        
    }
    
    #STORE results
    #my @ARRAY_info = split("_", $occ);
    #my $kindex = $ARRAY_info[1];
    if(!( -d "./kmasker_raw_plots_$PID")){
        mkdir("kmasker_raw_plots_$PID");
    }   
    my $LIST = new IO::File($list, "r") or die "\n unable to read $list $!";    
    while(<$LIST>){
        next if($_ =~ /^$/);
        next if($_ =~ /^#/);
        my $line = $_;
        $line =~ s/\n//;
        if(-e $line.".png"){
            system("mv ".$line.".png kmasker_raw_plots_$PID" . "/".$line."_hist.png");
        }
    }
    
    #clean    
    if(-e "kmasker_seq_$PID.ids"){
        system("rm kmasker_seq_$PID.ids");
    }
}

## subroutine
#	Direct generation of plots for user filterd input
#
sub plot_maker_direct{
	my $file		=	$_[0];
	my $type		=	$_[1];
	
	#HEADER
	my $systemcall_header 	= `head -n 1 $file`;
	my @ARRAY_header		= split("\t", $systemcall_header);
	if(scalar(@ARRAY_header) != 3){
		print "\n\n This does not look like the expected format (3 columns). Process stopped!\n\n";
		exit();
	}	
	
	#MAKE BOXPLOT
    if($type eq "boxplot"){
    	system("make_boxplot.R --out KMASKER_boxplot_".$loctime.".png --column1 ".$ARRAY_header[1]." --column2 ".$ARRAY_header[2]." -i ".$file);
    	system("make_boxplot.R --out KMASKER_boxplot_log_".$loctime.".png -l --column1 ".$ARRAY_header[1]." --column2 ".$ARRAY_header[2]." -i ".$file);
    }
    
    #MAKE heaxgon plot
    if($type eq "hexplot"){
    	system("Rscript --vanilla ".$path . "/make_hexplot.R ".$file);
    }
    
    #MAKE violin plot
    if($type eq "violin"){
    	system("Rscript --vanilla ".$path . "/make_violinplot.R ".$file);
    }    
}

## subroutine
#	Manages data processing for generation of plots
#
sub plot_maker{
    my $cfile		=	$_[0];
    my $type		=	$_[1];  
    my $href_list  	=	$_[2];
    
    if(defined $href_list){
    	#work with selection list
    	my %HASH_list = %{$href_list};
    	my $R_CFILE	 = new IO::File($cfile, "r") or die "\n unable to read list $cfile $!";
		my $W_CFILE	 = new IO::File("tmp_".$cfile, "w") or die "\n unable to write tmp_$cfile $!";			
		
		print "\n .. start processing ".$cfile;
		
		if(scalar(keys %HASH_list) > 0){
			print "\n .. working with provided ID list\n";		
			my $counter = 0;
			while(<$R_CFILE>){
				next if($_ =~ /^$/);
				next if($_ =~ /^#/);
				print "\n COUNTER : ".$counter++;
				my $line = $_;
				my $orig = $line;
				$line =~ s/\n//;
				$line =~ s/ /\t/g;
				my @ARRAY_file = split(/\t/, $line);
					
				#check file format
				if(scalar(@ARRAY_file) < 16){
					#
					print "\n .. WARNING in line (".scalar(@ARRAY_file)."): ".$orig."\n"; 
					print "\n This line was skipped because it does not look like expected format (by Kmasker).";
					print "\n LINE: ".$orig."!\n";
				}else{
					if(exists $HASH_list{$ARRAY_file[0]}){
						print $W_CFILE $orig;
					}	
				}				
			}
			#set new file as input
			$cfile = "tmp_".$cfile;
		}   	
    }
    
    #MAKE BOXPLOT
    if($type eq "boxplot"){
    	&plot_boxplot($cfile)
    }
    
    #MAKE heaxgon plot
    if($type eq "hexplot"){
    	print "\n\n\n Please prepare your data accordingly (data input should have 3 columns) and use --file for generating a hexplot!";
    	print "\n Process stopped!\n\n";
    }
    
    #MAKE violin plot
    if($type eq "violin"){
    	print "\n\n\n Please prepare your data accordingly (data input should have 3 columns) and use --file for generating a violin plot!";
    	print "\n Process stopped!\n\n";
    }
    
    
    #clean
    system("rm tmp_".$cfile) if(-e "tmp_".$cfile);
    
}


## subroutine
#
sub plot_violin{
    my $cfile	=	$_[0];
    my $list	=	$_[1];    
    
    #IMPLEMENT if costumer defined plot is required
    
}

## subroutine
#
sub plot_hexagon{
    my $cfile	=	$_[0];
    my $list	=	$_[1];
    
    #IMPLEMENT if costumer defined plot is required
        
}

## subroutine
#
sub plot_boxplot{
    my $cfile	=	$_[0];
       
    #edit header
    system("grep -v \"#\" ".$cfile." >C_".$cfile);
    my $W_CFILE	 = new IO::File("H_".$cfile, "w") or die "\n unable to write H_$cfile $!";			
    my @ARRAY_header = ("","","","","","","","","","","","","","","","");
    $ARRAY_header[0] =~ s/#//g;   	
  	$ARRAY_header[6] = "avg_kmer_count_in_set1";
    $ARRAY_header[7] = "avg_kmer_count_in_set2";
    $ARRAY_header[9] = "proportion_foldchange_in_set1";
    $ARRAY_header[11] = "proportion_foldchange_in_set2";
    $ARRAY_header[13] = "proportion_absent_in_set1";
    $ARRAY_header[15] = "proportion_absent_in_set2";
    print $W_CFILE join("\t", @ARRAY_header)."\n";
    system("cat H_".$cfile." C_".$cfile." >N_".$cfile);
    
    #average
    system("make_boxplot.R --out KMASKER_boxplot_avg_log_".$loctime.".png -l --column1 ".$ARRAY_header[6]." --column2 ".$ARRAY_header[7]." -i N_".$cfile);   
  	system("make_boxplot.R --out KMASKER_boxplot_avg_".$loctime.".png --column1 ".$ARRAY_header[6]." --column2 ".$ARRAY_header[7]." -i N_".$cfile);   
  	
  	#fold change
  	#system("make_boxplot.R --out KMASKER_boxplot_pfc_log_".$loctime.".png -l --column1 ".$ARRAY_header[9]." --column2 ".$ARRAY_header[11]." -i N_".$cfile);   
  	system("make_boxplot.R --out KMASKER_boxplot_pfc_".$loctime.".png --column1 ".$ARRAY_header[9]." --column2 ".$ARRAY_header[11]." -i N_".$cfile); 
  	
  	#absent
  	#system("make_boxplot.R --out KMASKER_boxplot_absent_log_".$loctime.".png -l --column1 ".$ARRAY_header[13]." --column2 ".$ARRAY_header[15]." -i N_".$cfile);   
  	system("make_boxplot.R --out KMASKER_boxplot_absent_".$loctime.".png --column1 ".$ARRAY_header[13]." --column2 ".$ARRAY_header[15]." -i N_".$cfile); 
  	
  	system("rm N_".$cfile." C_".$cfile." H_".$cfile)  
   
}


## subroutine
#  Routine is plooting the means of k-mer frequency per size-binned sequence
sub plot_barplot{
    my $file	=	$_[0];
    my $bin		= 	$_[1];
    if(defined $bin){
		system("make_barplot.R -c avg -b ".$bin." -i ".$file);   
    }else{
    	system("make_barplot.R -c avg -i ".$file);   
    } 
}
## subroutine
#
sub report_statistics{
	my $occ	=	$_[0];
	my $gff	=	$_[1];
	my $out = 	$_[2];
	
	my $outtag = "";
	if(defined $out){
		$outtag = "_${PID}_".$out;
	}

	#Statistics
    print "\n .. call statistic calculation\n"; #if(!defined $silent);
	my $path 		= dirname abs_path $0;
 	
 	system("$path/stats.R " . "--input " .$occ . " --gff " . $gff . " --class sequence" . " --out report_statistics_sequences". $outtag .".txt " . " >>$log 2>&1");
	print "\n finished calculating statistics for sequences \n";
	system("$path/stats.R "  . "--input " . $occ . " --gff " . $gff . " --class KRC" .    " --out report_statistics_KRC".$outtag.".txt " . " >>$log 2>&1");
    print "\n finished calculating statistics for KRCs \n";
         
    if((!(-e "report_statistics_sequences".$outtag.".txt")) || (!(-e  "report_statistics_KRC".$outtag.".txt" ))) {
    	print "\nSome statistics could not be calculated. The main reason for this is that there are no significant features in the gff file.\n";
    }
    else {
		system("$path/stats_overview.R " . " -s report_statistics_sequences".$outtag.".txt " . " -k report_statistics_KRC".$outtag.".txt >>$log 2>&1");
    }
}


## subroutine
#
sub custom_annotation{
	my $fasta		=	$_[0];
    my $gff			=	$_[1];
    my $feature     =   $_[2];
    my $href_DB 	=   $_[3]; 
    my $href_info	=	$_[4];
    my $href_path	=	$_[5];


   	my %HASH_path		= %{$href_path};
    my %HASH_info 	= %{$href_info};
    my %HASH_DB		= %{$href_DB};
    my $threads = $HASH_info{"threads"};
    my $temp_dir = $HASH_info{"temp_path"};
   	my $path_blast  = $HASH_path{"blastn"};
   	my $path_makeblastdb = $HASH_path{"makeblastdb"};

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
            system("$path_makeblastdb -in \"".$db."\" -dbtype nucl ");
            print("BLASTdb was built. You can use the path to your fasta just with -db in the future.\n");   
        }
    }
    else{
    	print "\n\n Parameter dbfasta or db have to be provided for custom annotation!";
    	print   "\n Kmasker was stopped.\n\n";
    	exit();
    }
     kmasker::filehandler::extract_feature_gff($fasta, $gff, $feature, $HASH_info{"temp_path"});
     if(exists $HASH_info{"user setting blast"}) {
        my $parameterstring = $HASH_info{"user setting blast"};
        system("$path_blast -db \"" . $db . "\" -query " . "${temp_dir}/selected_" . $fasta . " -num_threads ".$threads." -outfmt 6 " . $parameterstring . " -ungapped" . "  -out ${temp_dir}/kmasker_blast.txt");
     }
     else{
      	system("$path_blast -db \"" . $db . "\" -query " . "${temp_dir}/selected_" . $fasta . " -perc_identity 80 -evalue 0.1 -num_threads ".$threads." -outfmt 6 -ungapped" . " -out ${temp_dir}/kmasker_blast.txt");
     }
    kmasker::filehandler::add_annotation_to_gff($gff, "${temp_dir}/kmasker_blast.txt");
    #kmasker::functions::add_annotation($fasta, $db, $gff, $feature ,$href_info);
    (my $name,my $path,my $suffix) = fileparse($gff, qr/\.[^.]*/);
    if(-e "$path${name}_with_annotation${suffix}") {
    	print("\nIntegration of annotation information in GFF finished!\n");
      	unlink($gff);
      	system("mv" . " " . "$path/${name}_with_annotation${suffix}" . " " . $gff);
      	
      	my $substring = "_with_annotation_with_annotation";
      	if($gff =~ /$substring/){
      		my $gff_new = $gff;
      		$gff_new =~ s/_with_annotation_with_annotationd/_with_annotation/;
      		system("mv ".$gff." ".$gff_new);
      	}
    }
    else{
        print("An annotated GFF was not created ($path${name}_with_annotation${suffix}). Something went wrong!");
    }
	
}





1;
