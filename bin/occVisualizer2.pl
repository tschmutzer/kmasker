#!/usr/bin/perl -w

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;

#global variables
my $occ;

my $result = GetOptions ("occ=s"   				=> \$occ,		# provide the occ file
						 ); 	



my $occ_file = Bio::SeqIO->new( -format => "qual",
				-file => "$occ" );

# open R process
open R, "|/opt/Bio/R/3.3.0/bin/R --vanilla --silent > /dev/null";

# import ggplot
print R 'require(ggplot2);'; 
print R 'require(zoo);'; 

while( my $occ_obj = $occ_file->next_seq ) {

    my $id = $occ_obj->display_id;
    my $seq = $occ_obj->seq;

    if( $id eq "#FREAK" ) { next; }

    print "Plotting contig: " . $id . "\n";
    print R 'png("' . $id. '.png", width=2048, height=1024);';
    print R "occs <- c(" . join(',', split(' ',$seq )) . ');';
    print R "mean10 <- rollmean(occs, 500, align= 'center', fill=0);";
    print R "positions <- 1:length(occs);";
    print R "dataframe <- data.frame(positions, occs, mean10);";
    print R "colnames(dataframe) <- c('pos', 'occ', 'mean10');";

    my $command = "ggplot(data = dataframe, aes(x = pos, y = occ))";
    $command .= " + geom_line( aes(colour = occ), size=1.1)";
    $command .= " + scale_colour_gradient(low='blue', high='red')";
    $command .= " + geom_area( aes( x = pos, y = mean10), fill='blue')";
    $command .= " + geom_line( aes( x = pos, y = mean10), size=1.2)";
    $command .= " + labs(title = 'k-mer distribution of $id', x='positions (bp)', y='k-mer frequencies')";
    print R $command . ";";
    
# print to R device
    print R "dev.off();";
    print "DONE\n";
}
