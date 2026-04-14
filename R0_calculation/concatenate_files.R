#!/usr/bin/env Rscript
rm(list=ls())
library(tidyverse)

f_final = lapply(list.files(path = "tabs/", pattern = "*.txt", full.names = T),
function(f_name) { read.table(f_name, header = T, sep = "\t")})
f_final = do.call("bind_rows", f_final)
write.csv2(f_final, "r0_all.txt", row.names = F)