#!/usr/bin/env perl -w
use Getopt::Long;
use IO::File;
use POSIX; 


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# date:         2018_08_29
# last update:	2018_08_30
my $version = "0.0.1 rc180830";
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

my $occ;
my $rept=10;
my $outname_bed;
my $outname_info;
my $foldchange = 5;
my $frame;
my $foldchange_usr;
my $header_normal;
my $help;

GetOptions(	'occ=s' 	=>	\$occ,
			'out=s' 	=>	\$outname_bed,
			'rept=i'	=>	\$rept_usr,
			'help' 		=> 	\$help);

#######################
##
##	MAIN
##

if(!defined $occ){
	exit();	
}

if(defined $rept_usr){
	$rept = $rept_usr;
}


print "\n .. start building BED file";
&build_bed();
print "\n .. finished process\n\n";

## END MAIN
#######################

#######################
##
##	subroutine
##
sub build_bed(){
	
	if(!defined $outname_bed){
		my $suffix = $occ;
		my $infix = "_RT".$rept.".bed";
		$suffix =~ s/\.occ$/$infix/;
		$outname_bed  = $suffix;
	
	}
	
	$outname_info = $outname_bed;
	if($outname_bed =~ /bed$/){
		$outname_info =~ s/bed$/info/;
	}else{
		$outname_info .= ".info";
	}
	
	print "\n\n BED outfile = ".$outname_bed."\n";
	
	my $OCC 		= new IO::File($occ, "r") or die "\n unable to read $occ $!";	
	my $BEDFILE		= new IO::File($outname_bed, "w") or die "\n unable to write $outname_bed $!";
	my $BEDFILE_DIR	= new IO::File("direct_".$outname_bed, "w") or die "\n unable to write $outname_bed $!";
	my $HANDLE_info	= new IO::File($outname_info, "w") or die "\n unable to write $outname_info $!";
	
	my $cid			= "";
	my $seq_counter = 0;
	my $hgh_counter	= 0;
	my $low_counter = 0;
	my $number		= 1;
	my $lid			= "";
	my @ARRAY_occ 	= ();
	my @ARRAY_status= ();	#contains list of status		#SAME LENGTH
	my @ARRAY_HGH	= ();	#contains list of HGH LENGTH	#SAME LENGTH
	my @ARRAY_LOW	= ();	#contains list of LOW LENGTH	#SAME LENGTH
	my $running_status = "";
	my $regionlength= "";
	
	
	while(<$OCC>){
		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		$line = $_;
		$line =~ s/\n//;	
		$line =~ s/ $//g;
		
		if($line =~ m/^>/){
			#OCC HEADER
			if($number != 1){
				#last section
				if($running_status eq "H"){
					push(@ARRAY_HGH, $regionlength);
		#			print "\n LAST PUSH ".$regionlength." to HGH";
				}else{
					push(@ARRAY_LOW, $regionlength);
		#			print "\n LAST PUSH ".$regionlength." to LOW";
				}				
				#write
				&write_info($HANDLE_info, $cid, \@ARRAY_status, \@ARRAY_HGH, \@ARRAY_LOW, $seq_counter, $hgh_counter, $low_counter);
				&write_bed($BEDFILE, $cid, \@ARRAY_status, \@ARRAY_HGH, \@ARRAY_LOW);
				#Clean
				@ARRAY_status	= ();	
				@ARRAY_HGH		= ();	
				@ARRAY_LOW		= ();	
				$running_status = "";
				$regionlength	= "";	
				$seq_counter	= 0;	
				$hgh_counter	= 0;
				$low_counter	= 0;		
			}
			my @ARRAY = split(" ", $line);
			$cid = $ARRAY[0];
			$cid =~ s/^>//;
			$number++;
			
		}else{
			#OCC CONTENT
			my @ARRAY = split(" ", $line);	
			for(my $i=0;$i<scalar(@ARRAY);$i++){				
				
				my $status = "";	
				$seq_counter++;			
				#check if high or low
				if($ARRAY[$i]>=$rept){
					#high
					$status = "H";
					$hgh_counter++;
				}else{
					#low
					$status = "L";
					$low_counter++
				}
				
				if($running_status eq ""){
					#FIRST base
					$running_status = $status;
					$regionlength 	= 1;
					push(@ARRAY_status, $status);
				}else{
					#ALL OTHERs
					if($running_status eq $status){
						#SAME satus --> extend
						$regionlength++;
					}else{
						#DIFF satus --> SAVE and start new
						if($running_status eq "H"){
							push(@ARRAY_HGH, $regionlength);
						}else{
							push(@ARRAY_LOW, $regionlength);
						}
						
						#FRESH
						push(@ARRAY_status, $status);
						$regionlength 	= 1;
						$running_status = $status;
					}
				}				
			}						
		}#END of sequence (IF ELSE)
	}
	#last
	if($running_status eq "H"){
		push(@ARRAY_HGH, $regionlength);
	}else{
		push(@ARRAY_LOW, $regionlength);
	}
	&write_info($HANDLE_info, $cid, \@ARRAY_status, \@ARRAY_HGH, \@ARRAY_LOW, $seq_counter, $hgh_counter, $low_counter);
	&write_bed($BEDFILE, $cid, \@ARRAY_status, \@ARRAY_HGH, \@ARRAY_LOW);
	$HANDLE_info->close();
}


#######################
##
##	subroutine: Merging adjcacent repeat sections
##
sub down_size_merge(){
	
	#pct
	
	#to be done
	
}

#######################
##
##	subroutine
##
sub write_info(){
	my $handle 		= $_[0];
	my $cid_this 	= $_[1];
	my $aref_status = $_[2];
	my $aref_hgh 	= $_[3];
	my $aref_low 	= $_[4];
	my $count_seq 	= $_[5];
	my $count_hgh 	= $_[6];
	my $count_low 	= $_[7];
	
	my @AS = @{$aref_status};
	my @AH = @{$aref_hgh};
	my @AL = @{$aref_low};	
	
	print $handle $cid_this."\tlength=".$count_seq."; high=".$count_hgh."; low=".$count_low."; pct_repeats=".sprintf("%.2f", ($count_hgh*100/$count_seq))."\n";
	print $handle "+\n";
	print $handle "S\t".join(" ", @AS)."\n";
	print $handle "+\n";
	print $handle "H\t".join(" ", @AH)."\n";
	print $handle "+\n";
	print $handle "L\t".join(" ", @AL)."\n";
	print $handle "+\n";
	
	my $length_AS = scalar(@AS);
	my $length_AH = scalar(@AH);
	my $length_AL = scalar(@AL);
	
	#never have same length
	if($length_AH==$length_AL){
		if($length_AS==$length_AL){
			print "n .. somwething went wrong for (".$cid_this."). Working with ambigous data !\n\n";
		}
	}else{
		
	}
}


#######################
##
##	subroutine
##
sub write_bed(){
	my $handle 		= $_[0];
	my $cid_this 	= $_[1];
	my $aref_status = $_[2];
	my $aref_hgh 	= $_[3];
	my $aref_low 	= $_[4];
	
	my @AS = @{$aref_status};
	my @AH = @{$aref_hgh};
	my @AL = @{$aref_low};	
	
	my $start_status;
	if($AS[0] eq "H"){
		$start_status = "H";
	}else{
		$start_status = "L";
	}
	
	#NOTE : BED format: position in BED column 2 is 0-based ; position in BED column 3 is 1-based
	
	my $first = "";
	my $BED_start;
	my $BED_end;
	for(my $i=0;$i<@AH;$i++){	
		
		if($first eq ""){
		#FIRST	
			$first = 0;
			if($start_status eq "H"){
				$BED_start 	= 0;
				$BED_end	= $AH[$i];
				print $handle $cid_this."\t".$BED_start."\t".$BED_end."\n";
				$BED_start 	= $AH[$i]+$AL[$i];
			}else{
				$BED_start 	= $AL[$i];
				$BED_end	= ($BED_start + $AH[$i]);
				print $handle $cid_this."\t".$BED_start."\t".$BED_end."\n";
				$BED_start 	= $AH[$i]+$AL[$i];
			}
		}else{
		#REST	
		
			if($start_status eq "H"){
			#HIGHSTART
				$BED_end	= $BED_start + 	$AH[$i];
				print $handle $cid_this."\t".$BED_start."\t".$BED_end."\n";
				$BED_start 	= $BED_end	 +	$AL[$i] if(defined $AL[$i]);
			}else{
			#LOWSTART	
				$BED_start 	+= $AL[$i];
				$BED_end	= $BED_start + $AH[$i];
				print $handle $cid_this."\t".$BED_start."\t".$BED_end."\n";
				$BED_start 	= $BED_end;				
			}			
		}
	}		
}


1;
