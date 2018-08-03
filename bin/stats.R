#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(getopt, zoo, scales, stringr)
#library('getopt')
#use getopt for input file and help description
spec=matrix(c(
  'input', 'i', 1, "character",
  'help', 'h', 0, "logical",
  'gff', 'g', 1, "character",
  'class', 'c', 1, "character",
  'out', 'o', 1, "character"
), byrow=TRUE, ncol=4)
opt=getopt(spec)
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

feature<-opt$class
outname=paste("report_statistics_", feature, "_", opt$gff , ".tab", sep="" );
if ( !is.null(opt$out) ) {
  outname=opt$out;
}

my.read.lines=function(fname) {
  s = file.info( fname )$size 
  buf = readChar( fname, s, useBytes=T)
  strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
  #https://www.r-bloggers.com/faster-files-in-r/
}

options(scipen = 999)
#https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html
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
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}


####
#choose the feature (contig, KRC; KRR)
gff<-gffRead(opt$gff)
gff_feature<-gff[which(gff[,"feature"]==feature),]
#pharse quality-data
qual_file <- my.read.lines(file.path(opt$input))
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
    if(pos != opentags[length(opentags)]) {
      next_pos_in_text <- opentags[which(opentags == pos) + 1]
      #Generate a vector of numbers out of the characters
      #We select a substring out of qual_file between two lines starting with >
      #Then we split the string into a list with a whitespace as sperator
      #After that we turn the list into a vector
      #We convert the characters into numbers finaly.
      #Our OCC format specification definies 80 characters per line at maximum (but mostly they are about 30 long)
      #This check will be tricked out by a file which is not in the specification
      occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text[1] - 1)], " ")))
      occsub<-data.frame(Name = id, ID=getAttributeField(gff_feature[which(gff_feature[,"seqname"]==id),"attributes"], "ID") ,Start=gff_feature[which(gff_feature[,"seqname"]==id),"start"] , End=gff_feature[which(gff_feature[,"seqname"]==id),"end"], Length=abs(gff_feature[which(gff_feature[,"seqname"]==id),"end"] - gff_feature[which(gff_feature[,"seqname"]==id),"start"]-1) ,Min = 0, Max=0, Avg=0, "zero_pos"=0, "Q25"=0, "Q50"=0, "Q75"=0)
      for(i in 1:nrow(occsub)){
        occs_temp<-occs[(occsub[i,"Start"]):(occsub[i,"End"])]
        occsub[i,"Min"]=min(occs_temp)
        occsub[i,"Max"]=max(occs_temp)
        occsub[i,"Avg"]=mean(occs_temp)
        occsub[i,"zero_pos"]=length(occs_temp[occs_temp==0])
        occsub[i,"Q25"]=quantile(occs_temp, 0.25)[[1]]
        occsub[i,"Q50"]=quantile(occs_temp, 0.5)[[1]]
        occsub[i,"Q75"]=quantile(occs_temp, 0.75)[[1]]
      }
    }
    else{
      next_pos_in_text <- length(qual_file)
      occs<-as.numeric(unlist(strsplit(qual_file[(pos+1):(next_pos_in_text)], " ")))      
      occsub<-data.frame(Name = id, ID=getAttributeField(gff_feature[which(gff_feature[,"seqname"]==id),"attributes"], "ID") ,Start=gff_feature[which(gff_feature[,"seqname"]==id),"start"] , End=gff_feature[which(gff_feature[,"seqname"]==id),"end"], Length=abs(gff_feature[which(gff_feature[,"seqname"]==id),"end"] - gff_feature[which(gff_feature[,"seqname"]==id),"start"]-1) ,Min = 0, Max=0, Avg=0, "zero_pos"=0, "Q25"=0, "Q50"=0, "Q75"=0)
      for(i in 1:nrow(occsub)){
        occs_temp<-occs[(occsub[i,"Start"]):(occsub[i,"End"])]
        occsub[i,"Min"]=min(occs_temp)
        occsub[i,"Max"]=max(occs_temp)
        occsub[i,"Avg"]=mean(occs_temp)
        occsub[i,"zero_pos"]=length(occs_temp[occs_temp==0])
        occsub[i,"Q25"]=quantile(occs_temp, 0.25)[[1]]
        occsub[i,"Q50"]=quantile(occs_temp, 0.5)[[1]]
        occsub[i,"Q75"]=quantile(occs_temp, 0.75)[[1]]
      }
    }
    if(first == 0){
#	write.table(occsub, file=paste("report_statistics_", feature, "_", opt$gff , ".tab", sep="" ), quote=FALSE, row.names = FALSE, append=FALSE)
	write.table(occsub, file=outname, quote=FALSE, row.names = FALSE, append=FALSE)
	first<-1
    }
    else{
#	write.table(occsub, file=paste("report_statistics_", feature, "_", opt$gff , ".tab", sep="" ), quote=FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
	write.table(occsub, file=outname, quote=FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
    }
  }
}
