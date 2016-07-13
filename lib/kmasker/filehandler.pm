package kmasker::filehandler;
use Exporter qw(import);
use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_sequence read_occ);
our @EXPORT_OK = qw(read_sequence read_occ);

#all credit goes to http://code.izzid.com/2011/10/31/How-to-read-a-fasta-file-in-perl.html

sub read_sequence { #works for fasta 
   my ($fh, $seq_info) = @_;

   $seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0; 
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
         my $h = $_;    
         $h =~ s/^>//;  
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;   
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         s/\s+//g;  # remove any white space
         $seq_info->{seq} .= $_;
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;   
   }    
}


sub read_occ { #works for occ
   my ($fh, $seq_info) = @_;

   $seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0; 
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
         my $h = $_;    
         $h =~ s/^>//;  
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;   
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         #s/\s+//g;  # remove any white space # We do not want this, because we are using split with whitespaces
         $seq_info->{seq} .= $_;
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;   
   }    
}

sub sequence_length {
   $seqfile = $_[0]; #just the filename of the sequence file
   open(my $seqfh, "<", "$seqfile") or die "Can not open $file\n";
   (my $name,my $path,my $suffix) = fileparse($seqfile);
   my %seqdata;
   while (read_sequence($fasta, \%seqdata)) {
      print $seqdata{header} . " : " . length($seqdata{seq}) . "\n";
   }
}

1;