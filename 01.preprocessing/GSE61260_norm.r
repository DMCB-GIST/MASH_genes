rm(list = ls())
#
dir_work = "c:/Dropbox/ddasic/02.tues/08.TRPC6_hepatic_fibrosis/data/"
dir_geo = "c:/Dropbox/ddasic/geo/"
setwd(dir_work)
library(data.table);library(matrixStats);library(limma);library(WGCNA)
library(org.Hs.eg.db)
library(oligo)


#
id_expr = "GSE61260"
file = list.files(path = dir_geo, pattern = id_expr)

#
idx1 = grep("series", file)
sample_storage = paste0(dir_geo, file[idx1])
sample = fread(sample_storage, fill = T)
sample = data.frame(sample)
sample = data.frame(t(sample[c(33, 42:45), -1]))
colnames(sample) = c("id", "sex", "age", "bmi", "group")
storage_sample = paste0(id_expr, "_sample.txt")
write.table(sample, storage_sample, sep = "\t", row.names = F)

#
for(i in 1:length(sample$id)){
  print(i)
  id_i = sample$id[i]
  file = list.files(path = "c:/geo/GSE61260_RAW/", pattern = id_i)
  file.i = paste0("c:/geo/GSE61260_RAW/", file)
  expr.i = read.celfiles(file.i)
  expr.i = rma(expr.i)
  expr.i = exprs(expr.i)
  colnames(expr.i) = id_i
  if(i == 1){
    expr = expr.i
  } else {
    expr = cbind(expr, expr.i[rownames(expr), ,drop = F])
  }
}
write.table(expr, "temp.txt", sep = "\t")
expr = fread("temp.txt")
expr = data.frame(expr)

#
source("c:/01.research/bi/GPL11532_symbol.R")

#selecting a probe with maxmean value
rownames(expr) = paste0("s", (1:nrow(expr)))
expr_summary = collapseRows(expr[, -1], rowID = rownames(expr), rowGroup = expr[, 1])
expr_summary = expr_summary$datETcollapsed
expr_storage = paste0(id_expr, "_norm.txt")
write.table(expr_summary, expr_storage, sep = "\t")

