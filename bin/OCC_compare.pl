#!/usr/bin/perl -w
use Getopt::Long;
use IO::File;
use POSIX; 

my $version = "0.0.1";


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author:       Thomas Schmutzer
# date:         2016_08_10
# last update:	2016_11_30
# institute:    @IPK Gatersleben
# version:		0.0.02
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

my @ARRAY_occ;
my %HASH_index;
my $occ_file1;
my $occ_file2;
my $foldchange = 5;
my $frame;
my $foldchange_usr;
my $help;

GetOptions(	'occ1=s' 	=>	\$occ_file1,
			'occ2=s' 	=>	\$occ_file2,
			'fc=s'		=>	\$foldchange_usr,
			'frame=i'	=>	\$frame,
			'help' 		=> 	\$help);

#######################
##
##	MAIN
##

if((!defined $occ_file1)||(!defined $occ_file2)){	
	exit();	
}

if(defined $foldchange_usr){
	$foldchange = $foldchange_usr;
}

my $loctime = localtime;
$loctime = strftime('%Y%m%d_%H%M',localtime); ## outputs 1208171008

&read_OCC1();
&read_OCC2();
## END MAIN
#######################

#######################
##
##	subroutine
##
sub read_OCC1(){

	my $OCCin 	 = new IO::File($occ_file1, "r") or die "\n unable to read list $occ_file1 $!";	
	my $occline = "";
	my $oid 	= "";
	my $index	= 0;
	while(<$OCCin>){
		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		$line = $_;
		$line =~ s/\n//;	
		if($line =~ m/^>/){
			$line =~ s/^>//;
			if($oid ne ""){
				$HASH_index{$oid} = $index;
				$ARRAY_occ[$index]= $occline;
				$occline = "";
				my @linearray = split(" ", $line);
				$oid = $linearray[0];
				$index++;				
			}else{
				my @linearray = split(" ", $line);
				$oid = $linearray[0];
			}
		}else{
			if($occline eq ""){
				$occline .= $line;
			}else{
				$occline .= " ".$line;
			}
		}	
	}
	#last
	$HASH_index{$oid} = $index;
	$ARRAY_occ[$index]= $occline;
	$occline = "";
	my @linearray = split(" ", $line);
	$oid = $linearray[0];
	
}


#######################
##
##	subroutine
##
sub read_OCC2(){	
	
	my $OUT;
	if(defined $frame){
		$OUT = new IO::File("KMASKER_comparitive_analysis_FC".$foldchange."_frame".$frame."_".$loctime.".stats", "w") or die "\n unable to read list $!";	
	}else{
		$OUT = new IO::File("KMASKER_comparitive_analysis_FC".$foldchange."_".$loctime.".stats", "w") or die "\n unable to read list $!";	
	} 
	
	print $OUT "sid\tlength\tstart\tend\tsum1\tsum2\tavg1\tavg2\t#pos. with ".$foldchange."-fold increase D1\tpositions with ".$foldchange."-fold in D1 [%]";
	print $OUT "\t".$foldchange."-fold increase D2\tpositions with ".$foldchange."-fold in D2 [%]";
	print $OUT "\tpositions absent D1\tpositions absent D1 [%]\tpositions absent D2\tpositions absent D2 [%]\n";
	
	
	my $OCCin2 	 = new IO::File($occ_file2, "r") or die "\n unable to read list $occ_file2 $!";	
	my $occline = "";
	my $oid 	= "";
	while(<$OCCin2>){
		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		$line = $_;
		$line =~ s/\n//;	
		if($line =~ m/^>/){
			$line =~ s/^>//;
			if($oid ne ""){
				my $index_OCC1 		= $HASH_index{$oid};
				my $occline_OCC1 	= $ARRAY_occ[$index_OCC1];
				
				if(defined $frame){
					&compare_frames($oid, $occline_OCC1, $occline, $OUT);
				}else{
					&compare($oid, $occline_OCC1, $occline, $OUT);
				}
			}
			$occline 			= "";
			my @linearray 		= split(" ", $line);
			$oid 				= $linearray[0];			
		}else{
			if($occline eq ""){
				$occline .= $line;
			}else{
				$occline .= " ".$line;
			}
		}	
	}
	#last
	my $index_OCC1 		= $HASH_index{$oid};
	my $occline_OCC1 	= $ARRAY_occ[$index_OCC1];
	if(defined $frame){
		&compare_frames($oid, $occline_OCC1, $occline, $OUT);
	}else{
		&compare($oid, $occline_OCC1, $occline, $OUT);
	}
}


#######################
##
##	subroutine
##
sub compare(){	
	my $this_oid 		= $_[0];
	my $this_occline1 	= $_[1];
	my $this_occline2 	= $_[2];
	my $HANDLER 		= $_[3];
	
	my @ARRAY_occ1 = split(" ", $this_occline1);
	my @ARRAY_occ2 = split(" ", $this_occline2);
	
	print "\n WARNING: unqequal length at seqeunce ".$this_oid."\n" if(scalar(@ARRAY_occ1) != scalar(@ARRAY_occ2));
		
	my $SUM1 = 0;
	my $SUM2 = 0;
	my $NUM_positions_with_foldchange_D1 = 0;
	my $NUM_positions_with_foldchange_D2 = 0;
	my $NUM_positions_with_absent_kmer_D1 = 0;
	my $NUM_positions_with_absent_kmer_D2 = 0;
	
	for(my $i = 0; $i < scalar(@ARRAY_occ1);$i++){
		my $value1 = $ARRAY_occ1[$i];
		my $value2 = $ARRAY_occ2[$i];
		$SUM1 += $value1;
		$SUM2 += $value2;
		if($value1 < 1){
			$value1 = 1;
			$NUM_positions_with_absent_kmer_D1++;
		}
		if($value2 < 1){
			$value2 = 1;
			$NUM_positions_with_absent_kmer_D2++;
		}
		
		if($value1 == $value2){
			#equal
		}elsif($value1 > $value2){			
			if($value1 > ($value2 * $foldchange)){
				$NUM_positions_with_foldchange_D1++;
			}			
		}elsif($value2 > $value1){			
			if($value2 > ($value1 * $foldchange)){
				$NUM_positions_with_foldchange_D2++;
			}			
		}		
	}
	
	my $AVG1 = sprintf("%.2f", ($SUM1/scalar(@ARRAY_occ1)));
	my $AVG2 = sprintf("%.2f", ($SUM2/scalar(@ARRAY_occ1)));
	my $ABS1 = sprintf("%.2f", ($NUM_positions_with_absent_kmer_D1*100/scalar(@ARRAY_occ1)));
	my $ABS2 = sprintf("%.2f", ($NUM_positions_with_absent_kmer_D2*100/scalar(@ARRAY_occ2)));
	
	print $HANDLER $this_oid."\t".scalar(@ARRAY_occ1)."\t1\t".scalar(@ARRAY_occ1)."\t".$SUM1."\t".$SUM2."\t".$AVG1."\t".$AVG2."\t".$NUM_positions_with_foldchange_D1."\t".sprintf("%.2f", ($NUM_positions_with_foldchange_D1*100/scalar(@ARRAY_occ1)));
	print $HANDLER "\t".$NUM_positions_with_foldchange_D2."\t".sprintf("%.2f", ($NUM_positions_with_foldchange_D2*100/scalar(@ARRAY_occ2)));
	print $HANDLER "\t".$NUM_positions_with_absent_kmer_D1."\t".$ABS1."\t".$NUM_positions_with_absent_kmer_D2."\t".$ABS2."\n";
	
}


#######################
##
##	subroutine
##
sub compare_frames(){	
	my $this_oid 		= $_[0];
	my $this_occline1 	= $_[1];
	my $this_occline2 	= $_[2];
	my $HANDLER 		= $_[3];
	
	my @ARRAY_occ1 = split(" ", $this_occline1);
	my @ARRAY_occ2 = split(" ", $this_occline2);
	
	print "\n WARNING: unqequal length at seqeunce ".$this_oid."\n" if(scalar(@ARRAY_occ1) != scalar(@ARRAY_occ2));
	
	my $frame_count = 0;
	my $SUM1 = 0;
	my $SUM2 = 0;
	my $NUM_positions_with_foldchange_D1 = 0;
	my $NUM_positions_with_foldchange_D2 = 0;
	my $NUM_positions_with_absent_kmer_D1 = 0;
	my $NUM_positions_with_absent_kmer_D2 = 0;
	my $start = 1;
		
	for(my $i = 0; $i < scalar(@ARRAY_occ1);$i++){
		
		if($i==0){
			#INIT sequence
			my $default = "-";
			print $HANDLER $this_oid."\t".scalar(@ARRAY_occ1)."\t1\t".scalar(@ARRAY_occ1)."\t".$default."\t".$default."\t".$default."\t".$default."\t".$default."\t".$default;
			print $HANDLER "\t".$default."\t".$default."\t".$default."\t".$default."\t".$default."\t".$default."\n";	
		}
		
		if($frame_count == $frame){
			#PRINT
			my $AVG1 = sprintf("%.2f", ($SUM1/$frame_count));
			my $AVG2 = sprintf("%.2f", ($SUM2/$frame_count));
			my $ABS1 = sprintf("%.2f", ($NUM_positions_with_absent_kmer_D1*100/$frame_count));
			my $ABS2 = sprintf("%.2f", ($NUM_positions_with_absent_kmer_D2*100/$frame_count));			
			print $HANDLER $this_oid."\t".$frame_count."\t".$start."\t".$i."\t".$SUM1."\t".$SUM2."\t".$AVG1."\t".$AVG2."\t".$NUM_positions_with_foldchange_D1."\t".sprintf("%.2f", ($NUM_positions_with_foldchange_D1*100/$frame_count));
			print $HANDLER "\t".$NUM_positions_with_foldchange_D2."\t".sprintf("%.2f", ($NUM_positions_with_foldchange_D2*100/$frame_count));
			print $HANDLER "\t".$NUM_positions_with_absent_kmer_D1."\t".$ABS1."\t".$NUM_positions_with_absent_kmer_D2."\t".$ABS2."\n";
			
			#RESET
			$SUM1 = 0;
			$SUM2 = 0;
			$NUM_positions_with_foldchange_D1 = 0;
			$NUM_positions_with_foldchange_D2 = 0;
			$NUM_positions_with_absent_kmer_D1 = 0;
			$NUM_positions_with_absent_kmer_D2 = 0;
			$frame_count = 0;
			$start = $i+1;
		}
		$frame_count++;
		
		my $value1 = $ARRAY_occ1[$i];
		my $value2 = $ARRAY_occ2[$i];
		$SUM1 += $value1;
		$SUM2 += $value2;
	
		if($value1 < 1){
			$value1 = 1;
			$NUM_positions_with_absent_kmer_D1++;
		}
		if($value2 < 1){
			$value2 = 1;
			$NUM_positions_with_absent_kmer_D2++;
		}
		
		if($value1 == $value2){
			#equal
		}elsif($value1 > $value2){			
			if($value1 > ($value2 * $foldchange)){
				$NUM_positions_with_foldchange_D1++;
			}			
		}elsif($value2 > $value1){			
			if($value2 > ($value1 * $foldchange)){
				$NUM_positions_with_foldchange_D2++;
			}			
		}			
	}
	#LAST
	my $AVG1 = sprintf("%.2f", ($SUM1/$frame_count));
	my $AVG2 = sprintf("%.2f", ($SUM2/$frame_count));
	my $ABS1 = sprintf("%.2f", ($NUM_positions_with_absent_kmer_D1*100/$frame_count));
	my $ABS2 = sprintf("%.2f", ($NUM_positions_with_absent_kmer_D2*100/$frame_count));		
	print $HANDLER $this_oid."\t".$frame_count."\t".$start."\t".($start+$frame_count-1)."\t".$SUM1."\t".$SUM2."\t".$AVG1."\t".$AVG2."\t".$NUM_positions_with_foldchange_D1."\t".sprintf("%.2f", ($NUM_positions_with_foldchange_D1*100/$frame_count));
	print $HANDLER "\t".$NUM_positions_with_foldchange_D2."\t".sprintf("%.2f", ($NUM_positions_with_foldchange_D2*100/$frame_count));
	print $HANDLER "\t".$NUM_positions_with_absent_kmer_D1."\t".$ABS1."\t".$NUM_positions_with_absent_kmer_D2."\t".$ABS2."\n";
		
	
	
}

1;