#!/usr/bin/env Rscript
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})
setRepositories(graphics = F, ind=c(1,2,5,7))
if (!require("getopt")) install.packages("getopt")
library(getopt) 
if (!require("reshape2")) install.packages("reshape2")
library(reshape2) 
if (!require("scales")) install.packages("scales")
library(scales) 
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2) 
#use getopt for input file and help description
spec=matrix(c(
  'input', 'i', 1, "character",
  'help', 'h', 0, "logical",
  'column', 'c', 1, "character",
  'binsize', 'b', 2, "numeric",
  'dynamic', 'd', 0, "logical"
  
), byrow=TRUE, ncol=4)
opt=getopt(spec)
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}
options(scipen = 999)

column<-opt$column
##
data<-read.delim(opt$input)
##
if(is.null(opt$dynamic)) {
length_levels<-cut(data[,"length"], c(0,500,1000,5000, max(data[,"length"])), dig.lab = 4, right = F, include.lowest = T, labels=c(">0", ">500", ">1000", ">5000"))
} else {pacman::p_load(Hmisc); length_levels<-cut2(data[,"length"], g=5, m=length(data[,"length"])/5)}

##
data_melt<-data.frame(length_levels=length_levels, value=data[,column])
##
if(! is.null(opt$binsize)){
  binwidth<-opt$binsize
} else {binwidth<-50}
##
plot<-ggplot(data_melt, aes(x=value)) + 
  geom_histogram(binwidth = binwidth, aes(fill=length_levels), boundary=0)  + xlim(c(0,quantile(data_melt[,"value"], 0.99))) +
  theme_gray(base_size = 25) +
  xlab(paste("Parameter:",column, sep=" ")) + ylab("Count") + guides(fill=guide_legend(title="Sequence Length"))
png(paste(tools::file_path_sans_ext(opt$input), "_", column ,"_boxplot.png", sep=""), width=2048, height=1024)
print(plot)
dev.off()