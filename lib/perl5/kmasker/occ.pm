package kmasker::occ;
use Exporter qw(import);
#use jellyfish;
use File::Basename;
use kmasker::filehandler;
use POSIX qw/ceil/;
use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(normalize_occ apply_occ merge_occ);
our @EXPORT_OK = qw(make_occ normalize_occ apply_occ merge_occ);

my $version_PM_occ 	= "0.0.3 rc190212";

sub normalize_occ{
	my $file = $_[0];
	my $depth = $_[1];
	my %seqdata;
	(my $name,my $path,my $suffix) = fileparse($file, qr/\.[^.]*/);
	open(my $occ, "<", "$file") or die "Can not open $file\n";
	open(my $occ_norm, ">", $path.$name."_normalized.occ") or die "Can not write to " . $path.$name."_normalized.occ\n";
	if ($depth == 1) {
		print "Nothing to do! Depth was 1. \n";
		return;
	}
	else {
		while (read_occ($occ, \%seqdata)) {
			print $occ_norm ">".$seqdata{header}."\n";
			my @values = split /\s+/ , $seqdata{seq};
			foreach(@values) {
				#print $_ . ": $depth = " . ceil($_/$depth) . "\n" ;
				$_ = ceil($_/$depth);

				}
			#print $occ_norm "@values\n";
			my $whitespace = 0;
			for (my $i = 0; $i < scalar(@values); $i++) {
				if($whitespace) {
					print $occ_norm " " . $values[$i];
				}
				else{
					print $occ_norm $values[$i];
					$whitespace = 1;
				}
				if((($i+1) % 25 == 0) && ($i+1 != scalar(@values))) {
					print $occ_norm "\n";
					$whitespace = 0;
				}
			}
			print $occ_norm "\n";
		}
	}
	close($occ);
	close($occ_norm);
}

sub multi_occ{
	my $threshold = $_[0];
	my $fold_change  = $_[1];
	my $occ1 = $_[2];
	my $occ2 = $_[3];
	my $prefix = $_[4];
	(my $name1,my $path1,my $suffix1) = fileparse($occ1, qr/\.[^.]*/);
	(my $name2,my $path2,my $suffix2) = fileparse($occ2, qr/\.[^.]*/);
	open(my $occ1_f, "<", "$occ1") or die "Can not open $occ1\n";
	open(my $occ2_f, "<", "$occ2") or die "Can not open $occ2\n";
	my %occ_data_1;
	my %occ_data_2;
		
	#Feedback
	print "\n .. start processing both OCC files" ;#if(!defined $silent);
	
	#FILE HANDLER for first and second output
	open(my $first, ">", $path1 . "/" . $prefix . $name1 . ".tab");
	open(my $second, ">", $path2 . "/" . $prefix . $name2 . ".tab");	
	
	while(read_occ($occ1_f, \%occ_data_1)) {
		read_occ($occ2_f, \%occ_data_2);
		my @out_values;	
		my @occ1_values = split /\s+/, $occ_data_1{seq};
		my @occ2_values = split /\s+/, $occ_data_2{seq};
		if($occ_data_1{header} ne $occ_data_2{header}) {
			print "Warning: Headers in occ files are different! " . $occ_data_1{header} . " != " . $occ_data_1{header}  . "\n";
		}
		if(scalar(@occ1_values) != scalar(@occ2_values)) {
			die "Can not merge! The files differ in length.\n " . scalar(@occ1_values) . "!=" . scalar(@occ2_values) . "\n";
		}
		for (my $i = 0; $i < scalar(@occ1_values); $i++) {
			if(($occ1_values[$i] >= $threshold) || ($occ2_values[$i] >= $threshold)) {
				if(($occ2_values[$i] - $occ1_values[$i]) == 0){
					$out_values[$i] = 0;
				}
				elsif($occ1_values[$i] == 0) {
					$out_values[$i]="inf";
				}
				else {
					$out_values[$i] = ($occ2_values[$i] - $occ1_values[$i])/$occ1_values[$i];
				}
			}
			else {
				$out_values[$i] = 0;
			}
		}
		my $last = 0; #0 - uncompareable or not significant, 1 - first occ, 2 - second occ
	
		# WRITE for first and second tab
		for(my $i = 0; $i < scalar(@out_values); $i++) {
			if($out_values[$i] >= $fold_change) {
				if ($last != 2) {
					if($i != 0 && $last==1) {
						print $first "\t" . $i-1 . "\n";
					}
					print $second $occ_data_2{header} . "\t" . $i . "\t";
					$last = 2;
				}
							#print $last . " ";
							#print $i . " " ;

			}
			elsif(($out_values[$i] < 0) && ($out_values[$i] <= -1*(1/$fold_change))) {
					if ($last != 1) {
						if($i != 0 && $last==2) {
							print $second "\t" . $i-1 . "\n";
						}
						print $first $occ_data_1{header} . "\t" . $i . "\t";
						$last = 1;
					}
							#print $last . " ";
							#print $i . " " ;
				}
			elsif(($out_values[$i] == 0) || (($out_values[$i] > 0) && ($out_values[$i] < $fold_change)) || (($out_values[$i] < 0) && ($out_values[$i] > -1*(1/$fold_change)))) {
				if($last == 1) {
					print $first $i-1 . "\n";
					$last = 0;
					#print "closing frst";
				}
				elsif($last == 2) {
					print $second $i-1 . "\n";
					$last = 0;
				}
			}
			elsif($i == scalar(@out_values) -1 ) {
				if($last == 1) {
					print $first $i-1 . "\n";
					$last = 0;
					#print "closing frst";
				}
				elsif($last == 2) {
					print $second $i-1 . "\n";
					$last = 0;
				}
			}
			else{
				print "Internal error while processing!\n";
			}
		}	
	}
}

sub merge_occ {
	my $first_occ = $_[0];
	my $second_occ = $_[1];
	my $out_occ = $_[2];
	my %f_occdata;
	my %s_occdata;
	open(my $focc_stream, "<", "${first_occ}") or die "Can not open $first_occ\n";
	open(my $socc_stream, "<", "${second_occ}") or die "Can not open $second_occ\n";
	open(my $merged_occ, ">", "${out_occ}");
	while(read_occ($focc_stream, \%f_occdata)) {
		read_occ($socc_stream, \%s_occdata);
		my @focc_values = split /\s+/, $f_occdata{seq};
		my @socc_values = split /\s+/, $s_occdata{seq};
		if($f_occdata{header} ne $s_occdata{header}) {
			print "Warning: Headers in occ and fasta are different! " . $s_occdata{header} . " != " . $f_occdata{header}  . "\n";
		}
		if(scalar(@focc_values) != scalar(@socc_values)) {
			die "Can not merge! The files differ in length.\n " . scalar(@socc_values) . "!=" . scalar(@focc_values) . "\n";
		}
		print $merged_occ $f_occdata{header} . "\n";
		my @out_values;
		for (my $i = 0; $i < scalar(@focc_values); $i++) {
			$out_values[$i] = $focc_values[$i] + $socc_values[$i];
		}
		my $whitespace = 0;
		for (my $i = 0; $i < scalar(@out_values); $i++) {
			if($whitespace) {
				print $merged_occ " " . $out_values[$i];
			}
			else{
				print $merged_occ $out_values[$i];
				$whitespace = 1;
			}
			if((($i+1) % 25 == 0) && ($i+1 != scalar(@out_values))) {
				print $merged_occ "\n";
				$whitespace = 0;
			}
		}
		print $merged_occ "\n";
	}
}

1;