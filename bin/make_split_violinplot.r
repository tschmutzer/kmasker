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
## my $command = "Rscript --vanilla make_split_violinplot.R kmer_data_AB_small.txt";
## system($command);
## ####################

#!/usr/bin/env Rscript
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2) 
if (!require("data.table")) install.packages("data.table")
library(data.table) 
if (!require("dplyr")) install.packages("dplyr")
library(dplyr) 

#initialize split violin plot function
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})
geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

args = commandArgs(trailingOnly=TRUE)

# Test. Are input requirements fulfilled.
if (length(args)<1) {
  stop("No input data provided!", call.=FALSE)
} else if (length(args)==1) {
  # default output file

  # RUN
  
}

date<-format(Sys.Date(),"%y%m%d")
input<-args[1]
dt<-fread(input,h=T)
A<-colnames(dt)[2]
B<-colnames(dt)[3]
colnames(dt)<-c("kmer","A","B")
output<-paste0("Kmasker_plots_split_violin_",A,"_",B,"_",date,".png")

tab1<-data.table(m=A,y=dt[,2],x="x")
tab2<-data.table(m=B,y=dt[,3],x="x")
colnames(tab1)<-c("m","y","x")
colnames(tab2)<-c("m","y","x")
mtab<-rbind(tab1,tab2)
colnames(mtab)<-c("legend","frequency","distribution")

# HEADER of columns will be used as axis names
t<-ggplot(mtab,aes(distribution,frequency,fill=legend))+geom_split_violin()+theme_bw()+labs(title=paste0("k-mer comparison between ",A," and ",B))+theme(plot.title=element_text(hjust=0.5))

# generate visualisation in PNG format
png(output, width = 1800, height = 1800, res = 300)
t
dev.off()