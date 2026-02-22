rm(list = ls())
library(data.table)
library(org.Hs.eg.db)
library(WGCNA)
setwd("c:/01.research/02.tues/08.TRPC6_hepatic_fibrosis/data/")

###
access_id = "GSE135251"
sample = list.files(path = "c:/geo/", pattern = access_id)
idx1 = grep("_series", sample)
sample_storage = paste0("c:/geo/", sample[idx1])
sample = fread(sample_storage, fill = T)
sample = data.frame(sample)

###
sample = sample[c(33, 41, 42, 44, 45), -1]
sample = data.frame(t(sample))
colnames(sample) = c("id", "nas", "fibrosis", 
                     "dx", "stage")

###
file_sample = paste0(access_id, "_sample.txt")
write.table(sample, file = file_sample, row.names = F)

###
raw_dir = "c:/GEO/GSE135251_RAW/" 
id = sample$id

#
for (i in 1:length(id)) {
  print(i)
  id.i = id[i]
  raw_list = list.files(raw_dir, id.i)
  file_name = paste0(raw_dir, "/", raw_list)
  expr.i = read.table(file_name)
  colnames(expr.i) = c("ENSMBL", id.i)
  if (i == 1) {
    expr = expr.i
  } else {
    idx1 = match(expr[, 1], expr.i[, 1])
    expr = cbind(expr, expr.i[, 2, drop = F])
  }
}

#
expr[, 1] = mapIds(org.Hs.eg.db, keys = expr[, 1], keytype = "ENSEMBL", 
                   column = "SYMBOL")
expr = expr[!is.na(expr[, 1]), ]

###
rownames(expr) = paste0("s", (1:nrow(expr)))
expr_sum = collapseRows(expr[, -1], rowID = rownames(expr), rowGroup = expr[, 1])
datExpr = expr_sum$datETcollapsed

###
end_storage = paste0(access_id, "_count.txt")
write.table(datExpr, end_storage, sep = "\t")
