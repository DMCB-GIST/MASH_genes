rm(list = ls())
setwd("c:/01.research/02.tues/08.TRPC6_hepatic_fibrosis/data/")
library(data.table);library(matrixStats);library(limma);library(WGCNA)
library(org.Hs.eg.db)

#
id_expr = "GSE48452"
file = list.files(path = "c:/geo/", pattern = id_expr)
file = paste0("c:/geo/", file)
idx1 = grep("series", file)

#
expr = fread(file[idx1], fill = T)
expr = data.frame(expr)
names(expr)

#
id = unlist(expr[74, -1])
disease = unlist(expr[34, -1])
age = unlist(expr[43, -1])
age = gsub("age: ", "", age)
sex = unlist(expr[42, -1])
bmi = unlist(expr[44, -1])
fat = unlist(expr[38, -1])
inflammation = unlist(expr[39, -1])
nas = unlist(expr[45, -1])
fibrosis = unlist(expr[46, -1])
leptin = unlist(expr[50, -1])
adiponectin = unlist(expr[51, -1])
tissue = unlist(expr[36, -1])

sample = data.frame(id = id, disease = disease, age = age, sex = sex, bmi = bmi, fat = fat,
                    inflammation = inflammation, nas = nas, fibrosis = fibrosis, leptin = leptin, 
                    adiponectin = adiponectin, tissue)
sample_storage = paste0(id_expr, "_sample.txt")
write.table(sample, sample_storage, sep = "\t", row.names = F)

#GPL11532
#
expr = expr[-c(1:73), ]
write.table(expr, "temp.txt", sep = "\t", row.names = F, col.names = F)

#
expr = fread("temp.txt")
expr = data.frame(expr)
boxplot(expr[sample(1:nrow(expr), 1000), -1])
# expr[, -1] = normalizeQuantiles(expr[, -1])

#
source("c:/01.research/bi/GPL11532.R")

#selecting a probe with maxmean value
rownames(expr) = paste0("s", (1:nrow(expr)))
expr_summary = collapseRows(expr[, -1], rowID = rownames(expr), rowGroup = expr[, 1])
expr_summary = expr_summary$datETcollapsed
expr_storage = paste0(id_expr, "_norm.txt")
write.table(expr_summary, expr_storage, sep = "\t")

