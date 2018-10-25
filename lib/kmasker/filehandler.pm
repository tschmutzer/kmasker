package kmasker::filehandler;
use Exporter qw(import);
use strict;
use warnings;
use File::Basename;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_sequence read_occ sequence_length occ_length extract_sequence_region add_annotation_to_gff);
our @EXPORT_OK = qw(read_sequence read_occ sequence_length occ_length extract_sequence_region add_annotation_to_gff);

## VERSION
my $version_PM_filehandler 	= "0.0.1 rc170324";

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
         my @splitline = split(/(\s|\t)/ ,$h); #ignore any comments
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $splitline[0];
            return $seq_info;   
         }              
         else { # first time through only
            $seq_info->{header} = $splitline[0];
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
         my @splitline = split(/(\s|\t)/ ,$h); #ignore any comments
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $splitline[0];
            return $seq_info;   
         }              
         else { # first time through only
            $seq_info->{header} = $splitline[0];
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

sub extract_sequence_region {
   my $seqfile = $_[0];
   my $list = $_[1];
   my $offset = 0;
   if(defined $_[3]) {
      $offset = $_[3];
   }
   open(my $seqfh, "<", "$seqfile") or die "Can not open $seqfile\n";
   (my $name,my $path,my $suffix) = fileparse($seqfile, qr/\.[^.]*/);
   open(my $listfh, "<", "$list") or die "Can not open $list\n";
   my @start;
   my @stop;
   my @ident;
   while(<$listfh>) {
      my @line = split(/\t/, $_);
      push(@ident, $line[0]);
      push(@start, ($line[1] - $offset));
      push(@stop , ($line[2] - $offset));
   }
   my %seqdata;
   open(my $seql, ">", "$path/selected_$name$suffix") or die "Can not write to $path/selected_$name$suffix\n";
   my $oldi = 0;
   while (read_sequence($seqfh, \%seqdata)) {
      for(my $i = $oldi; $i < scalar(@ident); $i++) {
         if($seqdata{header} eq $ident[$i]) {
            print $seql ">" . $seqdata{header} . "_" . $start[$i] . "_" . $stop[$i] . "\n";
            my $temp = substr($seqdata{seq}, $start[$i], (($stop[$i] - $start[$i]) + 1));
            print $seql $temp . "\n";
         }
         if($seqdata{header} ne $ident[$i]) {
            $oldi = $i;
            last;
         }
      }
   }
}

sub extract_feature_gff {
   my $seqfile = $_[0];
   my $list = $_[1];
   my $feature = $_[2];
   my $temp = $_[3];
   my $offset = 1;
   open(my $seqfh, "<", "$seqfile") or die "Can not open $seqfile\n";
   (my $name,my $path,my $suffix) = fileparse($seqfile, qr/\.[^.]*/);
   if (defined $temp) {
      $path=$temp;
   }
   open(my $listfh, "<", "$list") or die "Can not open $list\n";
   my @start;
   my @stop;
   my @ident;
   while(<$listfh>) {
      if($_ !~ m/^#/ ) {
         my @line = split(/\t/, $_);
         #print $line[2];
         if($line[2] eq $feature) {
               push(@ident, $line[0]);
               push(@start, ($line[3] - $offset));
               push(@stop , ($line[4] - $offset));
         }
      }
   }
   my %seqdata;
   open(my $seql, ">", "$path/selected_$name$suffix") or die "Can not write to $path/selected_$name$suffix\n";
   my $oldi = 0;
   while (read_sequence($seqfh, \%seqdata)) {
      for(my $i = $oldi; $i < scalar(@ident); $i++) {
         if($seqdata{header} eq $ident[$i]) {
            print $seql ">" . $seqdata{header} . "_" . $start[$i] . "_" . $stop[$i] . "\n";
            my $temp = substr($seqdata{seq}, $start[$i], (($stop[$i] - $start[$i]) + 1));
            print $seql $temp . "\n";
         }
         if($seqdata{header} ne $ident[$i]) {
            $oldi = $i;
            last;
         }
      }
   }
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
   open( my $inFASTA, "<", "$fasta") or die "Can not open $fasta\n";
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
  my $tab=$_[0];
  my $length=$_[1];
  my $minlength=$_[2];
  my $featurename=$_[3];

   open( my $inTAB, "<", "$tab") or die "Can not read $tab (tab_to_gff) \n";
   open( my $lengthFH, "<", "$length") or die "Can not read $length (tab_to_gff) \n";

   (my $name,my $path,my $suffix) = fileparse($tab, qr/\.[^.]*/);
   open(my $outGFF, ">", $path . "/" . $name . ".gff") or die "Can not write gff (tab_to_gff) \n";
   print $outGFF "##gff-version 3". "\n";
   while(<$lengthFH>) {
      chomp($_);
      my @line = split(/\t/, $_);
            my $source = "Kmasker";
            my $type = "sequence";
            my $score = "."; #evalue
            my $strand = ".";
            my $phase = ".";
            my $ident = $line[0];
            my $start = 1;
            my $end = $line[1];
            my $attributes =  "ID=$ident;Name=$ident";
            print $outGFF  $ident . "\t" . $source . "\t" . $type . "\t" . $start . "\t" . $end . "\t" . $score  . "\t" . $strand . "\t" . $phase . "\t" . $attributes. "\n";
         }
   
   if (defined $_[4]) {
       my $subfeature = $_[4];
       my $subfeaturename = $featurename;
       if (defined $_[4]) {
          $subfeaturename = $_[5];
      }
      #print "\nUsing $subfeature as reference for subfeature annotation!\n";
      open (my $subTAB, "<", "$subfeature") or die "Can not read subfeature (tab_to_gff) \n";
      my $insub =  0;
      my $c = 1;
      my $s = 1;
      my $tabline = <$inTAB>;
      my $lengthcheck=1;
      while (<$subTAB>) {
         my @subline = split(/\t/, $_);
         my @line = split(/\t/, $tabline);
         my $ident = $line[0];
         my $start = $line[1] + 1;
         my $end = $line[2] + 1;
         my $ident_s = $subline[0];
         my $start_s = $subline[1] + 1;
         my $end_s = $subline[2] + 1;
         if(($end - $start) + 1 < $minlength) {
            $lengthcheck = 0;
         } 
         else{$lengthcheck = 1;}
         if($ident_s eq $ident && $start_s == $start && $end_s == $end){ #no subfeature
            my $source = "Kmasker";
            my $type = $featurename;
            my $score = "."; #evalue
            my $strand = ".";
            my $phase = ".";
            my $attributes =  "ID=${type}_$c;Name=${type}_$c;Parent=$ident";
            if($lengthcheck==1){
                  print $outGFF  $ident . "\t" . $source . "\t" . $type . "\t" . $start . "\t" . $end . "\t" . $score  . "\t" . $strand . "\t" . $phase . "\t" . $attributes. "\n";
                  $c++;
             }
             $tabline = <$inTAB>;
         }
         elsif($ident_s eq $ident && $start_s == $start && $end_s != $end) { #start of subfeature 
            $insub = 1;
            my $source = "Kmasker";
            my $type = $featurename;
            my $score = "."; #evalue
            my $strand = ".";
            my $phase = ".";
            my $attributes =  "ID=${type}_$c;Name=${type}_$c;Parent=$ident";
            if($lengthcheck==1){
               print $outGFF $ident . "\t" . $source . "\t" . $type . "\t" . $start . "\t" . $end . "\t"  . $score  . "\t" . $strand . "\t" . $phase . "\t" . $attributes. "\n";
            }
            #print subfeature
             $type = $subfeaturename;
             $score = "."; #evalue
             $strand = ".";
             $phase = ".";
             $attributes =  "ID=${type}_${c}.${s};Name=${type}_${c}.${s};Parent=${featurename}_${c}";
             if($lengthcheck==1){
               print $outGFF $ident_s . "\t" . $source . "\t" . $type . "\t" . $start_s . "\t" . $end_s . "\t"  . $score . "\t" . $strand . "\t" . $phase . "\t" . $attributes. "\n";
               $s++;
             }
         }
         elsif($ident_s eq $ident && $start_s != $start && $end_s == $end){ #end of subfeature
            if($insub == 1){
               $insub = 0;
               my $source = "Kmasker";
               my $type = $subfeaturename;
               my $score = "."; #evalue
               my $strand = ".";
               my $phase = ".";
               my $attributes =  "ID=${type}_${c}.${s};Name=${type}_${c}.${s};Parent=${featurename}_${c}";
               if($lengthcheck==1){
                   print $outGFF $ident_s . "\t" . $source . "\t" . $type . "\t" . $start_s . "\t" . $end_s . "\t"  . $score  . "\t" . $strand . "\t" . $phase .  "\t" . $attributes. "\n";
                   $c++;
               }
               $tabline = <$inTAB>;
               $s=0;
            }
            else{
               die "Error: Found end of subfeature without start!";
            }
         }
         elsif($insub == 1 && $ident_s eq $ident && $start_s > $start && $end_s < $end) { # in subfeature
               my $source = "Kmasker";
               my $type = $subfeaturename;
               my $score = "."; #evalue
               my $strand = ".";
               my $phase = ".";
               my $attributes =  "ID=${type}_${c}.${s};Name=${type}_${c}.${s};Parent=${featurename}_${c}";
               if($lengthcheck==1){
                  print $outGFF $ident_s . "\t" . $source . "\t" . $type . "\t" . $start_s . "\t" . $end_s . "\t" . $score  . "\t" . $strand . "\t" . $phase ."\t" . $attributes . "\n";
                  $s++;
               }
         }
         else{
            print "Warning: Internal error!\n Can not identify $ident_s $start_s $end_s to $ident $start $end!\n";
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
         my $source = "Kmasker";
         my $type = $featurename;
         my $score = "."; #evalue
         my $strand = ".";
         my $phase = ".";
         my $attributes =  "ID=${type}_$c;Name=${type}_$c;Parent=$ident";
         if(($end - $start) +1 < $minlength) { 
            print $outGFF $ident . "\t" . $source . "\t" . $type . "\t" . $start . "\t" . $end . "\t" . $score  . "\t" . $strand . "\t" . $phase . "\t" . $attributes . "\n";
         }
         $c++;
      }
   }  
}

sub merge_tab_seeds{ #check chomping ! 
   my $seeds = $_[0];
   #$second_seeds = $_[1];
   my $percent_length = $_[1];
   my $min = $_[2];
   open(my $seed_f, "<", "$seeds") or die "Can not open $seeds\n";
   my ($name1, $path1, $suffix1) = fileparse("$seeds", qr/\.[^.]*/);
   open(my $out , ">", $path1 . "/" . $name1 . "_Regions_merged.tab");
   my @ident;
   my @start;
   my @end;
   #my @output;
   my $temp_start=0;
   my $temp_stop=0;
   my $length = 0;

   while(<$seed_f>) {
      my @line = split(/\t/, $_);
      push(@ident, $line[0]);
      push(@start, $line[1]);
      push(@end, $line[2]);
      $length++;
   }
   my $old_start = 0;
   my $old_end = 0;
   for(my $i=1; $i<$length; $i++) {
      if($ident[$old_end] ne $ident[$i]){
         #push(@output, $ident[$old_end] . "\t" . $start[$old_start] . "\t" . $end[$old_end]);
         #if($start[$old_end] != $end[$old_end]) {
            print $out $ident[$old_end] . "\t" . $start[$old_start] . "\t" . $end[$old_end];
         #}
         $old_end = $i;
         $old_start = $i;
      }
      elsif(((($start[$i] - $end[$old_end]) / ((($end[$old_end] - $start[$old_end]) + 1 )+ (($end[$i] - $start[$i]) + 1))) <= $percent_length/100) || (($start[$i] - $end [$old_end]) <= $min )) {
         $old_end = $i;
      }
      else{
         #if($start[$old_end] != $end[$old_end]) {
             #push(@output, $ident[$old_end] . "\t" . $start[$old_start] . "\t" . $end[$old_end]);
             print $out $ident[$old_end] . "\t" . $start[$old_start] . "\t" . $end[$old_end];
         #}
         $old_end = $i;
         $old_start = $i;
      }
      #print $i . "\n";
   }
   #seperate output for the last element
   if($length > 0) {
      #push(@output, $ident[$old_end] . "\t" . $start[$old_start] . "\t" . $end[$old_end]);
      print $out $ident[$old_end] . "\t" . $start[$old_start] . "\t" . $end[$old_end];
   }
   #foreach (@output) {
    #  print $out "$_";
   #}
   close($out);
   close($seed_f);
}

sub add_annotation_to_gff{
   my $gff = $_[0];
   my $blast = $_[1];
   my $verbose = "true";
   if (defined $_[2]) {
      $verbose = $_[2];
   }
   open(my $gfffh, "<", "$gff") or die "Can not open $gff\n";
   open(my $blastfh, "<", "$blast") or die "Can not open $blast\n";
   (my $name,my $path,my $suffix) = fileparse($gff, qr/\.[^.]*/);
   open(my $gff_out, ">", "$path/${name}_with_annotation${suffix}") or die "Can not write to $path/${name}_with_annotation${suffix} \n";     
   my %blastresults;
   while(<$blastfh>){
      my @line = split(/\t/, $_);
      if(! exists $blastresults{$line[0]}) {
         $blastresults{$line[0]} = \@line;
      }
      else{
         print "Please run blast with -max_targer_seqs=1. I will skip all entries except the first one\n";
      }
   }
   print $gff_out "##gff-version 3\n";
   while(<$gfffh>) {
      if ($_ !~ m/^#/) {
         chomp($_);
         my @line = split(/\t/, $_);
         my $ident = $line[0] . "_" .($line[3] - 1) . "_" . ($line[4] - 1);
         if ( exists $blastresults{$ident}) {
            my $aliascount = () = $line[8] =~ /Alias/g;
            my $extension = "";
            if($aliascount > 0) {
               $extension = $aliascount + 1;
            }
            my $refname = @{$blastresults{$ident}}[1];
            $line[8] = $line[8] . ";Alias" . $extension . "=$refname";
            if ($verbose eq "true") {
               my $evalue = @{$blastresults{$ident}}[10];
               $evalue =~ s/^\s+|\s+$//g;
               my $score = @{$blastresults{$ident}}[11];
               $score =~ s/^\s+|\s+$//g;
               my $pident = @{$blastresults{$ident}}[3];
               $pident =~ s/^\s+|\s+$//g;
               my $blength = @{$blastresults{$ident}}[4];
               $blength =~ s/^\s+|\s+$//g;
               $line[8] = $line[8] . ";blast_score". $extension ."=$score;blast_evalue". $extension ."=$evalue;blast_identity". $extension ."=$pident;blast_algn_length". $extension ."=$blength";
            }
         }
         for(my $i=0; $i<scalar(@line)-1; $i++) {
          print $gff_out $line[$i] . "\t";
        }
        print $gff_out $line[scalar(@line)-1] . "\n";
      }
   }
}

1;
