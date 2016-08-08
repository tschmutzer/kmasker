#!/usr/bin/env Rscript
library('getopt')
#use getopt for input file and help description
spec=matrix(c(
  'input', 'i', 1, "character",
  'help', 'h', 0, "logical"
), byrow=TRUE, ncol=4)
opt=getopt(spec)
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}
#load zoo
library(zoo)
require(ggplot2)
require(scales)


#pharse quality-data
qual_file <- readLines(file(file.path(opt$input), open="r"))
#get the positions with >
opentags<-which(grepl(">.*", qual_file) == TRUE)

for (pos in opentags){
  #get the ID out of the whole line 
  id<-substring(qual_file[pos], 2, nchar(qual_file[pos])-1)
  #extract the numbers between two opening tags 
  #there should not be any comments (or other stuff) in the file, just numbers and opening tags
  if(pos < opentags[length(opentags)]) {
    next_pos_in_text <- opentags[which(opentags == pos) + 1]
    #Generate a vector of numbers out of the characters
    #We select a substring out of qual_file between two lines starting with >
    #Then we split the string into a list with a whitespace as sperator
    #After that we turn the list into a vector
    #We convert the characters into numbers finaly.
    occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text[1] - 1)], " "))) }
  else{
    next_pos_in_text <- length(qual_file)
    occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text)], " ")))
  }
  #this is the same part as in occVisualizer_v140902.pl
  print(paste("Plotting contig:", id, sep=" "));
  #png(paste(id, ".png", sep=""), width=2048, height=1024)
  mean10 <- rollmean(occs, 500, align= 'center', fill=0)
  positions <- 1:length(occs);
  dataframe <- data.frame(positions, occs, mean10)
  colnames(dataframe) <- c('pos', 'occ', 'mean10')
  plot <- ggplot(data = dataframe, aes(x = pos, y = occ)) + 
    geom_line( aes(colour = occ), size=1.1) +
    scale_colour_gradient(low='blue', high='red') +
    geom_area( aes( x = pos, y = mean10), fill='blue')  +
    geom_line( aes( x = pos, y = mean10), size=1.2)  +
    labs(title = paste("k-mer distribution of", id, sep=" "), x="position (bp)", y="k-mer frequency") +
    theme_gray(base_size = 19) + 
    theme(legend.text=element_text(size=15)) +
    guides(fill = guide_legend(keywidth = 3, keyheight = 2))
  png(paste(id, ".png", sep=""), width=2048, height=1024)
  #in a loop we have to explicitly print the plot
  print(plot)
  dev.off()
}