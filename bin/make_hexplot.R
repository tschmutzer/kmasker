#!/usr/bin/env Rscript

# ==========================================================
#
#
#     Kmasker plants visualistions         
#
#
# ==========================================================

## first version: 180727
## last update  : 180814

## ####################
## Perl (internal call)
## my $command = "Rscript --vanilla make_hexplot.R kmer_data_AB_small.txt";
## system($command);
## ####################

local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})
setRepositories(graphics = F, ind=c(1,2,5,7))

if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2) 
if (!require("data.table")) install.packages("data.table")
library(data.table) 
args = commandArgs(trailingOnly=TRUE)

# Test. Are input requirements fulfilled.
if (length(args)<1) {
  stop("No input data provided!", call.=FALSE)
} 

date<-format(Sys.Date(),"%y%m%d")
input<-args[1]
dt<-fread(input,h=T)
A<-colnames(dt)[2]
B<-colnames(dt)[3]
colnames(dt)<-c("k-mer","A","B")

output<-paste0("KMASKER_plots_hexplot_",A,"_",B,"_",date,".png")

# HEADER of columns will be used as axis names
p<-ggplot(dt,aes(A,B))+geom_hex(bins=30)+theme_classic()+labs(title=paste0("comparative analysis"),x=paste0(A),y=paste0(B))+theme(plot.title=element_text(hjust=0.5))

# generate visualisation in PNG format
png(output, width = 1800, height = 1800, res = 300)
p+geom_abline(intercept=0,col="red")
dev.off()
