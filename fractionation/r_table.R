#!/usr/bin/env Rscript

#make a table from all csv files in a directory

library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(readr)

multmerge = function(mypath){filenames = list.files(path=mypath, full.names=TRUE)
datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}

full_data = multmerge("/path/to/directory")

write.csv(full_data, file = "output.csv")
