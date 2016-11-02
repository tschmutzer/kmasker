package kmasker::filehandler;
use Exporter qw(import);
use strict;
use warnings;
use File::Basename;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_sequence read_occ sequence_length occ_length);
our @EXPORT_OK = qw(read_sequence read_occ sequence_length occ_length);

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
         if($seq_info->{seq}) {
            $seq_info->{seq} .= " " .$_;
         }
         else {
            $seq_info->{seq} = $_;
         }
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
   my $seqfile = $_[0]; #just the filename of the sequence file
   open(my $seqfh, "<", "$seqfile") or die "Can not open $seqfile\n";
   (my $name,my $path,my $suffix) = fileparse($seqfile, qr/\.[^.]*/);
   open(my $seql, ">", "$path/$name$suffix.length") or die "Can not write to $path/$name$suffix.length \n";
   my %seqdata;
   while (read_sequence($seqfh, \%seqdata)) {
      #print $seqdata{header} . " : " . length($seqdata{seq}) . "\n";
      print $seql $seqdata{header} . "\t" . length($seqdata{seq}) . "\n";
   }
   close($seqfh);
   close($seql);
}

sub occ_length {
   my $seqfile = $_[0]; #just the filename of the sequence file
   open(my $seqfh, "<", "$seqfile") or die "Can not open $seqfile\n";
   (my $name,my $path,my $suffix) = fileparse($seqfile, qr/\.[^.]*/);
   open(my $seql, ">", "$path/$name$suffix.length") or die "Can not write to $path/$name$suffix.length \n";
   my %seqdata;
   while (read_occ($seqfh, \%seqdata)) {
      #print $seqdata{header} . " : " . length($seqdata{seq}) . "\n";
      my @occvalues = split /\s+/, $seqdata{seq};
      print $seql $seqdata{header} . "\t" . scalar(@occvalues) . "\n";
   }
   close($seqfh);
   close($seql);
}

sub fasta_to_tab {
   my $fasta = $_[0];
   my $prefix = $_[1];
   # Initiating Handler 
   open( my $inFASTA, "<", "$fasta");
   (my $name,my $path,my $suffix) = fileparse($fasta, qr/\.[^.]*/);
   open( my $TAB, ">", $path . "/$prefix" . $name . ".tab");
      my %seqdata;   
      while(read_sequence($inFASTA, \%seqdata)) {  
        my $id         = $seqdata{header};
        my @ARRAY_id   = split(" ", $id);
        $id            = $ARRAY_id[0];   
        my $seq        = $seqdata{seq};
        while ($seq =~ m/(X|x)+/gc) {
           print $TAB $id .  "\t"  . $-[0] . "\t" . ($+[0] - 1) . "\n"; #needs at least perl 5.6.0
      }
                  
   }
   close($inFASTA);
   close($TAB);   
                  
}  

sub tab_to_gff {
  my $tab=$_[1];
  my $featurename=$_[0];
   open( my $inTAB, "<", "$tab");
   (my $name,my $path,my $suffix) = fileparse($tab, qr/\.[^.]*/);
   open(my $outGFF, ">", $path . "/" . $name . ".gff");
   print $outGFF "##gff-version 3". "\n";
   if (defined $_[2]) {
      my $subfeature = $_[2];
      print "Using $subfeature as reference for subfeature annotation!\n";
      open (my $subTAB, "<", "$subfeature");
      my $insub =  0;
      my $c = 1;
      my $s = 1;
      my $tabline = <$inTAB>;
      while (<$subTAB>) {
         my @subline = split(/\t/, $_);
         my @line = split(/\t/, $tabline);
         my $ident = $line[0];
         my $start = $line[1] + 1;
         my $end = $line[2] + 1;
         my $ident_s = $subline[0];
         my $start_s = $subline[1] + 1;
         my $end_s = $subline[2] + 1;
         if($ident_s eq $ident && $start_s == $start && $end_s == $end){ #no subfeature
            my $source = "kmasker";
            my $type = $featurename;
            my $score = "."; #evalue
            my $strand = "?";
            my $phase = ".";
            my $attributes =  "ID=${type}_$c;Name=${type}_$c";
            print $outGFF  $ident . "\t" . $source . "\t" . $type . "\t" . $start . "\t" . $end . "\t" . $score  . "\t" . $strand . "\t" . $phase . "\t" . $attributes. "\n";
             $tabline = <$inTAB>;
             $c++;
         }
         elsif($ident_s eq $ident && $start_s == $start && $end_s != $end) { #start of subfeature 
            $insub = 1;
            my $source = "kmasker";
            my $type = "main_".$featurename;
            my $score = "."; #evalue
            my $strand = "?";
            my $phase = ".";
            my $attributes =  "ID=${type}_$c;Name=${type}_$c";
            print $outGFF $ident . "\t" . $source . "\t" . $type . "\t" . $start . "\t" . $end . "\t"  . $score  . "\t" . $strand . "\t" . $phase . "\t" . $attributes. "\n";
            #print subfeature
             $type = "sub_".$featurename;
             $score = "."; #evalue
             $strand = "?";
             $phase = ".";
            $attributes =  "ID=main_${featurename}_${c}_${type}_${s};Name=${type}_${s};Parent=main_${featurename}_${c}";
            print $outGFF $ident_s . "\t" . $source . "\t" . $type . "\t" . $start_s . "\t" . $end_s . "\t"  . $score . "\t" . $strand . "\t" . $phase . "\t" . $attributes. "\n";
            $s++;
         }
         elsif($ident_s eq $ident && $start_s != $start && $end_s == $end){ #end of subfeature
            if($insub == 1){
               $insub = 0;
               my $source = "kmasker";
               my $type = "sub_".$featurename;
               my $score = "."; #evalue
               my $strand = "?";
               my $phase = ".";
               my $attributes =  "ID=main_${featurename}_${c}_${type}_${s};Name=${type}_${s};Parent=main_${featurename}_${c}";
               print $outGFF $ident_s . "\t" . $source . "\t" . $type . "\t" . $start_s . "\t" . $end_s . "\t"  . $score  . "\t" . $strand . "\t" . $phase .  "\t" . $attributes. "\n";
               $tabline = <$inTAB>;
               $c++;
               $s=0;
            }
            else{
               die "Error: Found end of subfeature without start!";
            }
         }
         elsif($insub == 1 && $ident_s eq $ident && $start_s > $start && $end_s < $end) { # in subfeature
               my $source = "kmasker";
               my $type = "sub_".$featurename;
               my $score = "."; #evalue
               my $strand = "?";
               my $phase = ".";
               my $attributes =  "ID=main_${featurename}_${c}_${type}_${s};Name=${type}_${s};Parent=main_${featurename}_${c}";
               print $outGFF $ident_s . "\t" . $source . "\t" . $type . "\t" . $start_s . "\t" . $end_s . "\t" . $score  . "\t" . $strand . "\t" . $phase ."\t" . $attributes . "\n";
               $s++;
         }
         else{
            print "Warning: Internal error!";
            $tabline = <$inTAB>;
            $c++;
            $s=0;
         }
      }
   }
   else{
      my $c = 1;
      while(<$inTAB>) {
         my @line = split(/\t/, $_);
         my $ident = $line[0];
         my $start = $line[1] + 1;
         my $end = $line[2] + 1;
         my $source = "kmasker";
         my $type = $featurename;
         my $score = "."; #evalue
         my $strand = "?";
         my $phase = ".";
         my $attributes =  "ID=${type}_$c;Name=${type}_$c";
         print $outGFF $ident . "\t" . $source . "\t" . $type . "\t" . $start . "\t" . $end . "\t" . $score  . "\t" . $strand . "\t" . $phase . "\t" . $attributes . "\n";
         $c++;
      }
   }  
}

1;