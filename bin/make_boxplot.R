#!/usr/bin/env Rscript
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})
setRepositories(graphics = F, ind=c(1,2,5,7)
if (!require("getopt")) install.packages("getopt", )
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
  'column1', '1', 1, "character",
  'column2', '2', 1, "character",
  'log', 'l', 0, "logical",
  'out', 'o', 0, "character"
), byrow=TRUE, ncol=4)
opt=getopt(spec)
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}
options(scipen = 999)

data<-read.delim(opt$input)
data_melt<-melt(data[,c(opt$column1, opt$column2)])
data_melt[,"value"]<-as.numeric(data_melt[,"value"])
data_melt[,"variable"]<-as.factor(data_melt[,"variable"])

plot<-ggplot(data_melt, aes(x=variable, y=value)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)  +
  theme_gray(base_size = 25)

if(!is.null(opt$log)) {
  outname<-paste(tools::file_path_sans_ext(opt$input), "_log", ".png", sep="")
  if(!is.null(opt$out)) {
    outname<-(opt$out)
  }
  plot<-plot +   labs(title = paste("boxplot of",levels(data_melt[,"variable"])[1], levels(data_melt[,"variable"])[2] ,sep=" ", collapse=" "), x="variable", y="value [log10]")
  plot<-plot+scale_y_log10(labels=trans_format("log10", math_format(10^.x)))
  png(outname, width=2048, height=1024)
  
} else {
  outname<-paste(tools::file_path_sans_ext(opt$input), ".png", sep="")
  if(!is.null(opt$out)) {
    outname<-(opt$out)
  }
  plot<-plot + labs(title = paste("boxplot of",levels(data_melt[,"variable"])[1], levels(data_melt[,"variable"])[2] ,sep=" ", collapse=" "), x="variable")
  png(outname, width=2048, height=1024)
  
}

#in a loop we have to explicitly print the plot
print(plot)
dev.off()
