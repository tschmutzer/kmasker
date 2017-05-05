#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use IO::File;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author:       Thomas Schmutzer
# date:         2012_08_25
my $last_update=	"2012_09_01";
# institute:    @IPK Gatersleben
my $version = "0.0.1rc002";
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




my $occ_file;;
my $occ_dir;
my $help;
my %H_occ_stats;	#key is occ_file
my %Hfrequencies;	#Hfrequencies{value}-->frequ

GetOptions(	'dir=s' => 	\$occ_dir,
			'occ=s' =>	\$occ_file,
			'help' 	=> 	\$help);

$help = 1 if((!defined $occ_dir)&&(!defined $occ_file));

if(defined $help){
	print "\n\n  version: ".$version."\n\n";
	exit();
}


#####################################################
##
##  MAIN
##
#reading directory
my @all_occ_files = ();
if(!defined $occ_file){
	if(!defined $occ_dir){
		$occ_dir = "." ;
	}elsif("." eq $occ_dir){
		
	}elsif(1 == $occ_dir){
		$occ_dir = "." ;
	}	
	
	opendir (DIR, $occ_dir) or die $!;	
	while (my $file = readdir(DIR)){
		push(@all_occ_files, $file) if($file =~ /\.occ$/); 
	}
}else{
	push(@all_occ_files, $occ_file);
}

#############
# computation round
my %H_results;
foreach my $occ_file (@all_occ_files){
	my $href = &init_stat();
	&compute_frequency($occ_file);
	my $href_loaded = &compute_statistics($href, $occ_file);
	$H_results{$occ_file} = $href_loaded;
}
&write_statistics(\%H_results);




######################################################


#######################
##
##	subroutine
##
sub init_stat(){
	################
	#definition F50:
	#	1.)	sorting all frequency values by value
	#	2.) F50 is the biggest value in the sorted 50% set (F50 would be median ... but all other Fx are named and defined by this)
	#
	my %H_stat = ();	
	$H_stat{"# zero"} 	= 0;
	$H_stat{"% zero"} 	= 0;
	$H_stat{"P5"} 		= 0;
	$H_stat{"P10"}	 	= 0;
	$H_stat{"P20"}	 	= 0;
	$H_stat{"P25"}	 	= 0;
	$H_stat{"P50"}	 	= 0;
	$H_stat{"P75"} 		= 0;
	$H_stat{"avg"} 	= 0;
	$H_stat{"median"} 	= 0;	#equal F50
	$H_stat{"%median2x"}= 0;
	$H_stat{"%avg2x"} 	= 0;
	$H_stat{"max"} 		= 0;
	$H_stat{"total"} 	= 0;
	return \%H_stat;	
}


#######################
##
##	subroutine
##
sub compute_frequency(){	
	
	my $occ = shift @_;
	my $OCCin 	= new IO::File($occ, "r") or die "could not open trace : ".$occ.": $!\n";
	%Hfrequencies = ();

	while(<$OCCin>){
		my $line = $_;
		next if($line =~ /^$/);
		next if($line =~ /^#/);
		next if($line =~ /^>/);
		$line =~ s/\n//;
		my @linearray = split(" ", $line);
		foreach(@linearray){
			if(exists $Hfrequencies{$_}){
				$Hfrequencies{$_}++;
			}else{
				$Hfrequencies{$_} = 1;
			}
		}
	}		
}


#######################
##
##	subroutine
##
sub compute_statistics(){
	my $href = shift @_;
	my $file = shift @_;
	
	my %H_STAT = %{$href};	
	my $total_sum = 0;
	my $total_position = 0;
	my $median;
	my $position_over_2median 	= 0;
	my $position_over_avg 		= 0;
	my $max;	
	
	my $FPcurve = new IO::File("PecentilCurve_".$file.".txt", "w") or die "could not open trace  $!\n";	
	my @array_of_FPcurve;
	
	#first run through hash
	foreach my $occ_value1 (keys %Hfrequencies){
		$total_position += $Hfrequencies{$occ_value1};
		$total_sum += $Hfrequencies{$occ_value1} * $occ_value1;
	}
	
	#print "ToSu ".$total_sum;
	
	$H_STAT{"# zero"} 	= $Hfrequencies{0} if(exists $Hfrequencies{0});
	$H_STAT{"% zero"} 	= sprintf("%.1f", (($Hfrequencies{0} * 100 ) / $total_position)) if(exists $Hfrequencies{0});
	$H_STAT{"avg"} 		= sprintf("%.1f", ($total_sum/$total_position)); 		

	
	my $running_position = 0;
	my $last_position	 = 0;	
	foreach my $occ_value (sort {$a <=> $b} keys %Hfrequencies){
		
		#print "\n".$occ_value.'\t'.$Hfrequencies{$occ_value};
		if(defined $max){
			#first element
			$max = $occ_value if($occ_value > $max);		
		}else{
			$max = $occ_value;
		}
		$running_position += $Hfrequencies{$occ_value};
	
	
		#compute & update FPcurve		
		my $positions2update = sprintf("%.0f",(($running_position * 100) / $total_position));
		for(my $i = $last_position; $i < $positions2update; $i++){
			$array_of_FPcurve[$i] = $occ_value;
		}
		$last_position = $positions2update;
				
		
		#update 2median
		if(defined $median){
			if($occ_value > (2 * $median)){
				$position_over_2median += $Hfrequencies{$occ_value} ;
				print "2M(median = ".$median." and occ = ".$occ_value.") (".$file."):".$Hfrequencies{$occ_value}." --> (".$position_over_2median.")\n";
			}
		}
		
		
		#median
		if($H_STAT{"median"} == 0){
			if($running_position > ($total_position / 2)){
				$H_STAT{"median"} 		= $occ_value;
				$median 				= $occ_value;
			}
			
		}		
		
		
		#update 2avg
		if($occ_value > (2*$H_STAT{"avg"})){
			$position_over_avg += $Hfrequencies{$occ_value};
		}
		#print "\n pos over ".$position_over_avg." \t ".$occ_value;
		
		#F75
		if($H_STAT{"P75"} == 0){
			$H_STAT{"P75"} = $occ_value if($running_position > ( 3 * ($total_position / 4)));	
		}
		
		#F50
		if($H_STAT{"P50"} == 0){
			$H_STAT{"P50"} = $occ_value if($running_position > ($total_position / 2));	
			$median = $occ_value;
		}
		
		#F25
		if($H_STAT{"P25"} == 0){
			$H_STAT{"P25"} = $occ_value if($running_position > ($total_position / 4));	
		}
		
		#F20
		if($H_STAT{"P20"} == 0){
			$H_STAT{"P20"} = $occ_value if($running_position > ($total_position / 5));	
		}
		
		#F10
		if($H_STAT{"P10"} == 0){
			$H_STAT{"P10"} = $occ_value if($running_position > ($total_position / 10));	
		}
		
		#F5
		if($H_STAT{"P5"} == 0){
			$H_STAT{"P5"} = $occ_value if($running_position > ($total_position / 20));	
		}
	}	
	
	$H_STAT{"%median2x"} 	= sprintf("%.1f", (($position_over_2median * 100 ) / $total_position));
	$H_STAT{"total"} 		= $total_position;
	$H_STAT{"max"} 			= $max;
	$H_STAT{"%avg2x"}		= sprintf("%.1f", (($position_over_avg * 100 ) / $total_position));
	
	
	#print DISTRIBUTION:
	print $FPcurve "pctg\tocc";
	my $pos = 0;
	foreach(@array_of_FPcurve){
		print $FPcurve "\n".$pos."\t".$_;
		$pos++;
	}
	$FPcurve->close();
	
	return \%H_STAT;
	
}



#######################
##
##	subroutine
##
sub write_statistics(){
	
	my $href_results = shift @_;
	my %H= %{$href_results};
	#write file
	my $OCCstat = new IO::File("OVERALL_occstatistics.txt", "w") or die "could not open trace  $!\n";
	
	#HEADER
	print $OCCstat "#name";
	print $OCCstat "\ttotal";
	print $OCCstat "\t# zero";
	print $OCCstat "\t% zero";
	print $OCCstat "\tavg";
	print $OCCstat "\tmedian";		#equal F50
	print $OCCstat "\t%median2x";
	print $OCCstat "\t%avg2x";
	print $OCCstat "\tmax";
	print $OCCstat "\tP5";
	print $OCCstat "\tP10";
	print $OCCstat "\tP20";
	print $OCCstat "\tP25";
	print $OCCstat "\tP50";
	print $OCCstat "\tP75";
	print $OCCstat "\n";
	
	my %ORDER = ();
	foreach my $key (sort keys %H){
		my %H_stat = %{$H{$key}};
		print $OCCstat $key;
		print $OCCstat "\t".$H_stat{"total"};
		print $OCCstat "\t".$H_stat{"# zero"};
		print $OCCstat "\t".$H_stat{"% zero"};
		print $OCCstat "\t".$H_stat{"avg"};
		print $OCCstat "\t".$H_stat{"median"};		#equal F50
		print $OCCstat "\t".$H_stat{"%median2x"};
		print $OCCstat "\t".$H_stat{"%avg2x"};
		print $OCCstat "\t".$H_stat{"max"};
		print $OCCstat "\t".$H_stat{"P5"};
		print $OCCstat "\t".$H_stat{"P10"};
		print $OCCstat "\t".$H_stat{"P20"};
		print $OCCstat "\t".$H_stat{"P25"};
		print $OCCstat "\t".$H_stat{"P50"};
		print $OCCstat "\t".$H_stat{"P75"};
		print $OCCstat "\n";		
	}	
	
	$OCCstat->close();
	
}

