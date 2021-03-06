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
if (!require("scales")){install.packages("scales")}
library("scales")
if (!require("stringr")){install.packages("stringr")}
library("stringr")

#use getopt for input file and help description
spec=matrix(c(
  'sequence_stats', 's', 1, "character",
  'help', 'h', 0, "logical",
  'KRC_stats', 'k', 1, "character",
  'pid', 'p', 2, "character"
), byrow=TRUE, ncol=4)
opt=getopt(spec)
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}
options(scipen = 999)
stats<-read.table(opt$sequence_stats, header=TRUE)
stats_krc<-read.table(opt$KRC_stats, header=TRUE)

###overview stats
###############
ncontigs<-nrow(stats)
nkrc<-nrow(stats_krc)
avgk10<-length(which(stats[,"Avg"]>10))
avgsequencelength<-mean(stats[,"Length"])
sumsequencelength<-sum(stats[,"Length"])
avgsequenceq50<-mean(stats[,"Q50"])
avgkrcq50<-mean(stats_krc[,"Q50"])
avgk<-mean(stats[,"Avg"])
avgkkrc<-mean(stats_krc[,"Avg"])
################


fileConn<-file(paste("KMASKER_report_overview_statistics_", opt$pid, ".txt", sep=""))
writeLines(c(
paste("Input file SEQ :", "\t", opt$sequence_stats ),
paste("Input file KRC :", "\t", opt$KRC_stats),
paste("Number of sequences:", "\t", ncontigs),
paste("Number of Kmer Repeat Clusters:", "\t", nkrc),
paste("Number of sequences with average Kmer density > 10:", "\t", avgk10),
paste("Total sequence length", "\t", sumsequencelength), 
paste("Average sequence length:", "\t", avgsequencelength), 
paste("Average Q50 of Kmer density in sequences :", "\t", avgsequenceq50), 
paste("Average Kmer density in sequences:", "\t", avgk), 
paste("Average Kmer density in Kmer Repeat Clusters:", "\t", avgkkrc)), fileConn)
close(fileConn)

