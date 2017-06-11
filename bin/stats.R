#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(getopt, zoo, scales, stringr)
#library('getopt')
#use getopt for input file and help description
spec=matrix(c(
  'input', 'i', 1, "character",
  'help', 'h', 0, "logical",
  'gff', 'g', 1, "character",
  'class', 'c', 1, "character" 
), byrow=TRUE, ncol=4)
opt=getopt(spec)
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}
options(scipen = 999)
gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

####
#choose the feature (contig, KRC; KRR)
feature<-opt$class
gff<-gffRead(opt$gff)
gff_feature<-gff[which(gff[,"feature"]==feature),]
#pharse quality-data
qual_file <- readLines(file(file.path(opt$input), open="r"))
#get the positions with >
opentags<-which(grepl(">.*", qual_file) == TRUE)
#if( !is.null(opt$list)) {
#  id_list<-scan(file = opt$list, what=character())
#} else {
#  id_list<-unlist(lapply(str_split(str_trim(str_sub(qual_file[opentags],2)), " "), `[[`, 1))
#}
#id_list<-str_trim(id_list)
id_list<-gff_feature[,"seqname"]
first<-0
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
     # approx<-(next_pos_in_text-pos)*30
     # if(approx > 300000) {
    #    print("Sorry, your contig seems to be too large to be read in. Anyway, you can force me to read it with --force")
    #    next;
   #   }
      #Our OCC format specification definies 80 characters per line at maximum (but mostly they are about 30 long)
      #This check will be tricked out by a file which is not in the specification
      occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text[1] - 1)], " ")))
    }
    else{
      next_pos_in_text <- length(qual_file)
    #  approx<-(next_pos_in_text-pos)*30
     # if(approx > 300000) {
    #    print("Sorry, your contig seems to be too large to be read in. Anyway, you can force me to read it with --force")
     #   next;
   #   }
      occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text)], " ")))
    }
    output<-data.frame(Name = id , Start=gff_feature$start, End=gff_feature$end, Min = min(occs), Max=max(occs), Avg=mean(occs), "zero_pos"=length(occs[occs==0]), "Q25"=quantile(occs, 0.25)[[1]], "Q50"=quantile(occs, 0.25)[[1]], "Q75"=quantile(occs, 0.25)[[1]])
    #write output here
    print(output)
    if(first == 0){
      write.table(output, file=paste("Stats_", feature, "_", opt$gff , ".tab", sep="" ), quote=FALSE, row.names = FALSE, append=FALSE)
      first<-1
    }
    else{
      write.table(output, file=paste("Stats_", feature, "_", opt$gff , ".tab", sep="" ), quote=FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
    }
  }
}