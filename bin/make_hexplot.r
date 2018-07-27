# ==========================================================
#
#
#     Kmasker plants visualistions         
#
#
# ==========================================================

## 180727

## ####################
## Perl (internal call)
## my $command = "Rscript --vanilla make_hexplot.R kmer_data_AB_small.txt";
## system($command);
## ####################

#!/usr/bin/env Rscript
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
colnames(dt)<-c("k-mer","data set A"," data set B")

output<-paste0("Kmasker_plots_hexplot_",A,"_",B,"_",date,".png")

# HEADER of columns will be used as axis names
p<-ggplot(dt,aes(A,B))+geom_hex(bins=30)+theme_classic()+labs(title=paste0("k-mer comparison between ",A," & ",B),x=paste0("k-mer count ",A),y=paste0("k-mer count ",B))+theme(plot.title=element_text(hjust=0.5))

# generate visualisation in PNG format
png(output, width = 1800, height = 1800, res = 300)
p
dev.off()