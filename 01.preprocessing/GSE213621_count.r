rm(list = ls())
library(data.table)
library(org.Mm.eg.db)
library(Orthology.eg.db)
library(org.Hs.eg.db)
library(WGCNA)
setwd("c:/Dropbox/ddasic/02.tues/08.TRPC6_hepatic_fibrosis/data/")

#
id_expr = "GSE213621"
files = list.files(path = "c:/Dropbox/ddasic/geo/", pattern = id_expr)

#
idx1 = grep("series", files)
sample = fread(paste0("c:/Dropbox/ddasic/geo/", files[idx1]), fill = T)
sample = data.frame(sample)
sample = data.frame(t(sample[c(29, 38, 39), -1]))
colnames(sample) = c("id", "tissue", "fibrosis")
sample_storage = paste0(id_expr, "_sample.txt")
write.table(sample, sample_storage, sep = "\t", row.names = F)

#
idx1 = grep("FPKMs", files)
expr = fread(paste0("c:/Dropbox/ddasic/geo/", files[idx1]))
expr = data.frame(expr)
idx1 = grep("TRP", expr$V1)

#
source("c:/Dropbox/ddasic/bi/gene_summarization.r")
gene_num = table(expr[, 1])
gene_con = names(gene_num[gene_num == 1])
gene_dup = names(gene_num[gene_num != 1])

#
idx1 = match(gene_con, expr[, 1])
expr_con = expr[idx1, ]
expr_con = data.frame(expr_con, row.names = 1)

#
idx1 = which(expr[, 1] %in% gene_dup)
expr = expr[idx1, ]

#selecting a probe with maxmean value
rownames(expr) = paste0("s", (1:nrow(expr)))
expr_summary = collapseRows(expr[, -1], rowID = rownames(expr), rowGroup = expr[, 1])
expr_summary = expr_summary$datETcollapsed
expr_summary = rbind(expr_con, expr_summary)
#expr_storage = paste0(id_expr, "_norm.txt")
expr_storage = paste0(id_expr, "_count.txt")
write.table(expr_summary, expr_storage, sep = "\t")
