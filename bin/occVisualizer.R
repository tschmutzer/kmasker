#!/usr/bin/env Rscript
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})
setRepositories(graphics = F, ind=c(1,2,5,7))
if (!require("getopt")){install.packages("getopt")}
library("getopt")
if (!require("zoo")){install.packages("zoo")}
library("zoo")
if (!require("pillar")){install.packages("pillar")}
library("pillar")
if (!require("ggplot2")){install.packages("ggplot2")}
library("ggplot2")
if (!require("scales")){install.packages("scales")}
library("scales")
if (!require("stringr")){install.packages("stringr")}
library("stringr")
#use getopt for input file and help description
spec=matrix(c(
  'input', 'i', 1, "character",
  'help', 'h', 0, "logical",
  'force', 'f', 0, "logical",
  'sws', 'w', 2, "numeric",
  'list', 'l', 2, "character",
  'log', 'g', 0, "logical",
  'dynamic', 'd', 0, "logical"
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
wsize=500
if( !is.null(opt$sws)) {
  wsize=opt$sws
  print(paste("Sliding window was set to ", wsize, sep=""))
}
if (!is.null(opt$dynamic)) {
  dynamic=2
  print("Dynamic sliding window estimation activated");
}

my.read.lines=function(fname) {
  s = file.info( fname )$size 
  buf = readChar( fname, s, useBytes=T)
  strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
  #https://www.r-bloggers.com/faster-files-in-r/
}

#pharse quality-data
qual_file <- my.read.lines(file.path(opt$input))
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
        print("Sorry, your sequence seems to be too large to be read in. Anyway, you can force me to read it with --force")
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
        print("Sorry, your sequence seems to be too large to be read in. Anyway, you can force me to read it with --force")
        next;
      }
      occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text)], " ")))
    }
    #this is the same part as in occVisualizer_v140902.pl
    print(paste("Plotting contig:", id, sep=" "));
    if(length(occs)>1000000){
      if(force == TRUE) {
        print("Warning: Your sequence has more than 1000000bp. The process will take some time!")
      } else {
        print(paste("Finally your sequence has", length(occs), "bp. This is too much to print. You can force me with --force", sep=" "))
        next;
      }
    }
    if (!is.null(opt$dynamic)) {
      temp<-length(occs)*dynamic/100
      exp<-trunc(log10(temp))
      temp<-temp/(10^exp)
      values<-c(1,2,5)
      if(length(which(temp-values>0))==0) {
        print("Sequence too short for dynamic window estimation.\n");next;
      }
      value<-max(values[which(temp-values>0)])
      wsize<-value*10^exp
    }
    #png(paste(id, ".png", sep=""), width=2048, height=1024)
    mean10 <- rollmean(occs, wsize, align= 'center', fill=0)
    m<-mean(occs)
    q25<-quantile(occs, 0.25)[[1]]
    q50<-quantile(occs, 0.5)[[1]]
    q75<-quantile(occs, 0.75)[[1]]
    q90<-quantile(occs, 0,9)[[1]]
    positions <- 1:length(occs);
    dataframe <- data.frame(positions, occs, mean10)
    dataframe<-as.data.frame(apply(dataframe, c(1,2), function(x) if(x<1){return(1)}else{return(x)}))
    colnames(dataframe) <- c('pos', 'occ', 'mean10')
     plot <- ggplot(data = dataframe, aes(x = pos, y = occ)) + 
        geom_area( aes( x = pos, y = mean10), fill='blue')  +
        geom_line( aes( x = pos, y = mean10), size=1.2) 
    if(! is.null(opt$log)){
        plot + labs(title = paste("k-mer distribution of", id, "sws", wsize, sep=" "), x="position (bp)", y="k-mer frequency [log10]")
        plot + scale_y_log10(breaks=trans_breaks("log10",function(x) 10^x), labels=trans_format("log10", math_format(10^.x))) +

    } else {
      plot <- ggplot(data = dataframe, aes(x = pos, y = occ)) + 
        plot + labs(title = paste("k-mer distribution of", id, "sws", wsize, sep=" "), x="position (bp)", y="k-mer frequency") +
    }
    plot + theme_gray(base_size = 19) + theme(legend.text=element_text(size=15))
    values<-c()
      if(m > 0) { 
        plot + geom_line(aes(y=m, colour="avg"), linetype=3, size=1.2)
        values<-c(values, "avg" = "blue")
      }
      if(q25 > 0) { 
        plot + geom_line(aes(y=q25, colour="q25"), linetype=3, size=1.2)
        values<-c(values, "q25" = "red")
      }
      if(q50 > 0) { 
        plot + geom_line(aes(y=q50, colour="q50"), linetype=3, size=1.2)
        values<-c(values, "q50" = "brown")
      }
      if(q75 > 0) { 
        plot + geom_line(aes(y=q75, colour="q75"), linetype=3, size=1.2)
        values<-c(values, "q75" = "darkgreen")
      }      
      if(m > 0) { 
        plot + geom_line(aes(y=q90, colour="q90"), linetype=3, size=1.2)
        values<-c(values, "q90" = "violet")
      }
      plot + scale_color_manual(values = values) +labs(color="")
    png(paste(id, ".png", sep=""), width=2048, height=1024)
    #in a loop we have to explicitly print the plot
    print(plot)
    dev.off()
  }
}
