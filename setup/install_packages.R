#!/usr/bin/env Rscript

local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})
setRepositories(graphics = F, ind=c(1,2,5,7))
if (!require("getopt")){install.packages("getopt")}
if (!require("zoo")){install.packages("zoo")}
if (!require("pillar")){install.packages("pillar")}
if (!require("ggplot2")){install.packages("ggplot2")}
if (!require("scales")){install.packages("scales")}
if (!require("stringr")){install.packages("stringr")}
if (!require("reshape2")){install.packages("reshape2")}
if (!require("data.table")){install.packages("data.table")}
if (!require("dyplr")){install.packages("dyplr")}

