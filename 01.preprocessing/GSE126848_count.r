rm(list = ls())
library(data.table)
library(org.Hs.eg.db)
library(WGCNA)
setwd("c:/Dropbox/ddasic/02.tues/08.TRPC6_hepatic_fibrosis/data/")

###
access_id = "GSE126848"
files = list.files(path = "c:/Dropbox/ddasic/geo/", pattern = access_id)

#
idx1 = grep("raw", files)
file_storage = paste0("c:/Dropbox/ddasic/geo/", files[idx1])

#
expr = fread(file_storage, header = T)
expr = data.frame(expr)

#
idx1 = grep("_series", files)
sample_storage = paste0("c:/Dropbox/ddasic/geo/", files[idx1])
sample = fread(sample_storage, fill = T)
sample = data.frame(sample)

###
sample = sample[c(43, 38), -1]
sample = data.frame(t(sample))
colnames(sample) = c("id", "dx")

###
file_sample = paste0(access_id, "_sample.txt")
write.table(sample, file = file_sample, row.names = F)

#
expr[, 1] = mapIds(org.Hs.eg.db, keys = expr[, 1], keytype = "ENSEMBL", 
                   column = "SYMBOL")
expr = expr[!is.na(expr[, 1]), ]

source("c:/Dropbox/ddasic/bi/gene_summarization.r")
id_expr = access_id
expr_storage = paste0(id_expr, "_count.txt")
write.table(expr_summary, expr_storage, sep = "\t")
