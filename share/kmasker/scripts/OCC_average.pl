#!/usr/bin/env perl -w
use Getopt::Long;
use IO::File;
use POSIX; 


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# date:         2018_08_13
# last update:	2018_08_29
my $version = "0.0.4 rc180829";
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

my $occ;
my $max;
my $outname;
my $foldchange = 5;
my $frame;
my $foldchange_usr;
my $header_normal;
my $help;

GetOptions(	'occ=s' 	=>	\$occ,
			'out=s' 	=>	\$outname,
			'max=s'		=>	\$max,
			'help' 		=> 	\$help);

#######################
##
##	MAIN
##

if(!defined $occ){
	exit();	
}


print "\n .. start calculating averages";
&calculate_avg();
print "\n .. finished calculation";

## END MAIN
#######################

#######################
##
##	subroutine
##
sub calculate_avg(){
	
	if(!defined $outname){
		my $name = $occ;
		$name =~ s/^KMASKER_//;
		$outname = "KMASKER_average_".$name;
		if(defined $max){
			$outname = "KMASKER_average_m".$max."_".$name;
		}
	}
	
	my $suffix = $outname;
	$suffix =~ s/occ$/txt/;
	$outname = $suffix;

	my $OCC 	 = new IO::File($occ, "r") or die "\n unable to read $occ $!";	
	my $OUT 	 = new IO::File($outname, "w") or die "\n unable to write $outname $!";	
	print $OUT "ID\tlength\tavg\n";
	my $counter = 0;
	my $summer  = 0;
	my $number	= 0;
	my $lid		= "";
	my @ARRAY_occ = ();
	
	while(<$OCC>){
		next if($_ =~ /^$/);
		next if($_ =~ /^#/);
		$line = $_;
		$line =~ s/\n//;	
		$line =~ s/ $//g;
		
		if($line =~ m/^>/){
			#OCC HEADER
			
			if($number != 0){
				#not first --> calcutale
				my $count = scalar(@ARRAY_occ);
				foreach (@ARRAY_occ){
 					$summer += $_;
				}
				my $avg = sprintf("%.0f", ($summer/$count));
				if(defined $max){
					print $OUT $lid."\t".$count."\t".$avg."\n" if($avg <= $max);
				}else{
					print $OUT $lid."\t".$count."\t".$avg."\n";
				}
				
				#reset
				@ARRAY_occ = ();	
				$summer = 0;
				$count = 0;
			}
			
			$line =~ s/^>//;
			my @ARRAY_tmp = split(" ", $line);
			$lid = $ARRAY_tmp[0];
			
			$number++;
		}else{
			#OCC CONTENT
			my @ARRAY_tmp = split(" ", $line);
			push(@ARRAY_occ, @ARRAY_tmp);			
		}
	}
	#last
	my $count = scalar(@ARRAY_occ);
	foreach (@ARRAY_occ){
 		$summer += $_;
	}
	my $avg = sprintf("%.0f", ($summer/$count));
	if(defined $max){
		print $OUT $lid."\t".$count."\t".$avg."\n" if($avg <= $max);
	}else{
		print $OUT $lid."\t".$count."\t".$avg."\n";
	}	
}

1;
