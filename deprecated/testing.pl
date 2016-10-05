use jellyfish;

my $qf = jellyfish::QueryMerFile->new(shift(@ARGV));
foreach my $m (@ARGV) {
  my $mer = jellyfish::MerDNA->new($m);
  #$mer->canonicalize;
  print($mer, " ", $qf->get($mer), "\n");
}