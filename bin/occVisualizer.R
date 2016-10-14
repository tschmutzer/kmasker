#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(getopt, zoo, ggplot2, scales, stringr)
#library('getopt')
#use getopt for input file and help description
spec=matrix(c(
  'input', 'i', 1, "character",
  'help', 'h', 0, "logical",
  'force', 'f', 0, "logical",
  'list', 'l', 2, "character"
), byrow=TRUE, ncol=4)
opt=getopt(spec)
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}
options(scipen = 999)
#load zoo
#library(zoo)
#require(ggplot2)
#require(scales)
force=FALSE
if( !is.null(opt$force)) {
  force=TRUE
}

#pharse quality-data
qual_file <- readLines(file(file.path(opt$input), open="r"))
#get the positions with >
opentags<-which(grepl(">.*", qual_file) == TRUE)
if( !is.null(opt$list)) {
  id_list<-scan(file = opt$list, what=character())
} else {
  id_list<-unlist(lapply(str_split(str_trim(str_sub(qual_file[opentags],2)), " "), `[[`, 1))
}
id_list<-str_trim(id_list)

for (pos in opentags){
  if(str_split(str_trim(str_sub(qual_file[pos],2)), " ")[[1]][1] %in% id_list) {
    #get the ID out of the whole line 
    id<-str_split(str_trim(str_sub(qual_file[pos],2)), " ")[[1]][1]
    #extract the numbers between two opening tags 
    #there should not be any comments (or other stuff) in the file, just numbers and opening tags
    if(pos < opentags[length(opentags)]) {
      next_pos_in_text <- opentags[which(opentags == pos) + 1]
      #Generate a vector of numbers out of the characters
      #We select a substring out of qual_file between two lines starting with >
      #Then we split the string into a list with a whitespace as sperator
      #After that we turn the list into a vector
      #We convert the characters into numbers finaly.
      approx<-(next_pos_in_text-pos)*30
      if(approx > 300000 & force == FALSE) {
        print("Sorry, your contig seems to be too large to be read in. Anyway, you can force me to read it with --force")
        next;
      }
      #Our OCC format specification definies 80 characters per line at maximum (but mostly they are about 30 long)
      #This check will be tricked out by a file which is not in the specification
      occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text[1] - 1)], " ")))
      }
    else{
      next_pos_in_text <- length(qual_file)
      approx<-(next_pos_in_text-pos)*30
      if(approx > 300000 & force == FALSE) {
        print("Sorry, your contig seems to be too large to be read in. Anyway, you can force me to read it with --force")
        next;
      }
      occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text)], " ")))
    }
    #this is the same part as in occVisualizer_v140902.pl
    print(paste("Plotting contig:", id, sep=" "));
    if(length(occs)>1000000){
      if(force == TRUE) {
        print("Warning: Your contig has more than 1000000bp. The process will take some time!")
      } else {
        print(paste("Finally your Contig has", length(occs), "bp. This is too much to print. You can force me with --force", sep=" "))
        next;
      }
    }
    #png(paste(id, ".png", sep=""), width=2048, height=1024)
    mean10 <- rollmean(occs, 500, align= 'center', fill=0)
    positions <- 1:length(occs);
    dataframe <- data.frame(positions, occs, mean10)
    dataframe<-as.data.frame(apply(dataframe, c(1,2), function(x) if(x<1){return(1)}else{return(x)}))
    colnames(dataframe) <- c('pos', 'occ', 'mean10')
    plot <- ggplot(data = dataframe, aes(x = pos, y = occ)) + 
      geom_line( aes(colour = occ), size=1.1) +
      scale_colour_gradient(low='blue', high='red', guide = guide_colorbar(title = "frequency\n")) +
      geom_area( aes( x = pos, y = mean10), fill='blue')  +
      geom_line( aes( x = pos, y = mean10), size=1.2)  +
      labs(title = paste("k-mer distribution of", id, sep=" "), x="position (bp)", y="k-mer frequency [log10]") +
      theme_gray(base_size = 19) + 
      theme(legend.text=element_text(size=15)) +
      guides(fill = guide_legend(keywidth = 3, keyheight = 2)) +
      scale_y_log10(breaks=trans_breaks("log10",function(x) 10^x), labels=trans_format("log10", math_format(10^.x)))
    png(paste(id, ".png", sep=""), width=2048, height=1024)
    #in a loop we have to explicitly print the plot
    print(plot)
    dev.off()
  }
}