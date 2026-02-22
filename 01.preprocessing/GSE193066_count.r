rm(list = ls())
library(data.table)
library(org.Mm.eg.db)
library(Orthology.eg.db)
library(WGCNA)
setwd("c:/Dropbox/ddasic/02.tues/08.TRPC6_hepatic_fibrosis/data/")

#
id_expr = "GSE193066"
files = list.files(path = "c:/Dropbox/ddasic/geo/", pattern = id_expr)

#
idx1 = grep("series", files)
sample = fread(paste0("c:/Dropbox/ddasic/geo/", files[idx1]), fill = T)
sample = data.frame(t(sample[c(30, 41, 42, 43, 44), -1]))
colnames(sample) = c("id", "age", "sex", "fibrosis", "nafld_score")
sample$fibrosis = gsub("fibrosis stage: ", "", sample$fibrosis)
sample$fibrosis = as.numeric(sample$fibrosis)
storage_sample = paste0(id_expr, "_sample.txt")
write.table(sample, storage_sample, sep = "\t", row.names = F)

#
idx1 = grep("164", files)
expr = fread(paste0("c:/Dropbox/ddasic/geo/", files[idx1]), skip = 2)
expr = data.frame(expr)
idx1 = match(sample$id, colnames(expr))
idx1 = idx1[!is.na(idx1)]
expr = expr[, c(1, idx1)]

expr = data.frame(expr, row.names = 1)

expr_storage = paste0(id_expr, "_count.txt")
write.table(expr, expr_storage, sep = "\t")
