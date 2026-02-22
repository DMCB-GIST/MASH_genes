rm(list = ls())
library(data.table)
library(org.Mm.eg.db)
library(Orthology.eg.db)
library(WGCNA)
library(DESeq2)
library(gplots)
setwd("c:/Dropbox/ddasic/05.fri/01.NASH_gene/data/")

#
id_expr = "GSE48452"
files = list.files(pattern = id_expr)

#
idx1 = grep("norm", files)
expr = fread(files[idx1])
expr = data.frame(expr, row.names = 1)

#
#source("c:/Dropbox/ddasic/unsupervied_analysis/heatmap_genes.r")

#
idx1 = grep("sample", files)
sample = fread(files[idx1])
sample = data.frame(sample)

#
idx1 = match(colnames(expr), sample$id)
sample = sample[idx1, ]

#
phenotype = sample$disease
table(phenotype)
pheno_color = c("darkolivegreen1", "gold3", "darkred")
pheno_color = c("darkolivegreen1", "darkolivegreen3", "darkred", "grey")
pheno_type = c("CN", "Obese", "NASH", "Others")

#
idx.c = grep("Control", phenotype)
idx.c1 = grep("Healthy", phenotype)
idx.d1 = grep("Nash", phenotype)
idx.ro = grep(",", phenotype)
idx.c = setdiff(idx.c, idx.ro)
idx.c1 = setdiff(idx.c1, idx.ro)
idx.d1 = setdiff(idx.d1, idx.ro)

#
phenotype[] = pheno_type[4]
phenotype[idx.c] = pheno_type[1]
phenotype[idx.c1] = pheno_type[2]
phenotype[idx.d1] = pheno_type[3]
names(pheno_color) = pheno_type

# #
# idx1 = which(phenotype == "Others")
# expr = expr[, -idx1]
# sample = sample[-idx1]
# phenotype = phenotype[-idx1]

#
rowSds = rowSds(as.matrix(expr))
cut_off = quantile(rowSds, probs = 0.9)
idx1 = which(rowSds > cut_off)
expr = expr[idx1, ]
boxplot(expr)


expr_raw = expr
expr = t(expr)
expr = scale(expr)
expr = t(expr)
quantile(expr)
ext = 5
expr[expr >= ext] = ext
expr[expr <= -1*ext] = -1*ext

#
library(circlize)
library(ComplexHeatmap)
col = list(pheno_color)
names(col) = "Dx"
HA = HeatmapAnnotation(Dx = phenotype, 
                       col = col, 
                       which = "column",
                       simple_anno_size = unit(1, "cm"), 
                       show_legend = T)

#
col_heatmap = colorRamp2(c(-5, -1, 0, 1, 5), 
                         c("#428bca", "#6fcb9f", "grey95", "#ffe28a", "#fb2e01"))
figure_storage = paste0("F.heatmap_", id_expr, ".pdf")
pdf(figure_storage, height = 6*0.7, width = 10*0.7)

#
lgd = Legend(col_fun = col_heatmap, 
             title = "Z-score",
             direction = "horizontal")

#“topleft”, “topcenter”, “leftcenter”, “lefttop”, “leftcenter-rot”, “lefttop-rot”
ht <- Heatmap(expr,  
              show_row_dend = FALSE, 
              show_column_dend = FALSE, 
              col = col_heatmap,
              top_annotation = HA, 
              name = "Z-score", 
              show_row_names = FALSE, 
              show_column_names = FALSE, 
              show_heatmap_legend = F)

#
draw(ht, padding = unit(c(2, 2, 15, 2), "mm"))
draw(lgd, x = unit(0.01, "npc"), y = unit(0.99, "npc"), just = c("left", "top"))

dev.off()
