rm(list = ls())
setwd("c:/01.research/02.tues/08.TRPC6_hepatic_fibrosis/data/")
library(data.table);library(WGCNA);library(limma)

#
id_expr = "GSE159676"
file = list.files(path = "c:/geo/", pattern = id_expr)
file = paste0("c:/geo/", file)

#
idx1 = grep("series", file)
expr = fread(file[idx1], fill = T)
expr = data.frame(expr)
expr = expr[-nrow(expr), ]

#
sample = t(expr[c(37, 61), -1])
sample = data.frame(sample)
colnames(sample) = c("group", "id")

#
sample_storage = paste0(id_expr, "_sample.txt")
write.table(sample, sample_storage, sep = "\t", row.names = F)

#
expr = expr[-c(1:60), ]
write.table(expr, "temp.txt", sep = "\t", col.names = F, row.names = F)

#
expr = fread("temp.txt")
expr = data.frame(expr)
boxplot(expr[sample(1:nrow(expr), 1000), -1])
expr[, -1] = log(expr[, -1]+1)
expr[, -1] = normalizeQuantiles(expr[, -1]+1)

#
source("c:/01.research/bi/GPL6244.R")

#selecting a probe with maxmean value
rownames(expr) = paste0("s", (1:nrow(expr)))
expr_summary = collapseRows(expr[, -1], rowID = rownames(expr), rowGroup = expr[, 1])
expr_summary = expr_summary$datETcollapsed
expr_storage = paste0(id_expr, "_norm.txt")
write.table(expr_summary, expr_storage, sep = "\t")
