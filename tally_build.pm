package lib::tally_build;
use Exporter qw(import);
use File::Basename;
use strict;
use warnings;

our @ISA = qw(Exporter);
our @EXPORT = qw(
setup_index
make_occ_tally
make_kindex_tally
);
our @EXPORT_OK = qw(setup_index);


sub setup_index{
	my ($file, $k) = @_;
	system("head -n 100 ".$file." >input_head.txt");
	system("tail -n 100 ".$file." >input_tail.txt");
	system("cat input_head.txt input_tail.txt >input.txt");	
	my $md5_name = `md5sum input.txt`;
	system("rm input_head.txt input_tail.txt input.txt");
	
	#create index
	&make_kindex_tally($file, $md5_name, $k);
	
}


sub make_occ_tally{
	my ($fasta, $md5_name, $k) = @_;
	
	#
	&get_md5_name_from_repository();
	
	#3
	open PROC, "gt tallymer search -output qseqnum qpos counts -tyr ".$md5_name."_".$k." -q ".$fasta."|";
    my @result = <PROC>;
    chomp(@result);
	
}


sub make_kindex_tally{
	
	my ($input_seq, $md5_name, $k) = @_;
	
	#1
	#system("/opt/Bio/genometools/1.5.8/bin/gt suffixerator -dna -pl -tis -suf -lcp -db SMALL_ERR418039_2_trim_Q20.fastq.fasta -indexname SMALL");
	system("/opt/Bio/genometools/1.5.8/bin/gt suffixerator -dna -pl -tis -suf -lcp -db $input_seq -indexname ".$md5_name."_".$k);
	#2
	system("/opt/Bio/genometools/1.5.8/bin/gt tallymer mkindex -esa ".$md5_name."_".$k." -mersize 21 -minocc 1 -pl -indexname ".$md5_name."_".$k." -counts");
	
}

1;