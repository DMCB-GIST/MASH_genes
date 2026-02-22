rm(list = ls())
library(data.table)
library(org.Mm.eg.db)
library(Orthology.eg.db)
library(WGCNA)
library(DESeq2)
library(gplots)
setwd("c:/Dropbox/ddasic/05.fri/01.NASH_gene/data/")

#
id_expr = "GSE126848"
files = list.files(pattern = id_expr)

#
idx1 = grep("_count", files)
expr = fread(files[idx1])
expr = data.frame(expr, row.names = 1)
expr = log(expr + 1)

#
#source("c:/Dropbox/ddasic/unsupervied_analysis/heatmap_genes.r")

#
idx1 = grep("_sample", files)
sample = fread(files[idx1])
sample = data.frame(sample)

#
phenotype = sample$dx
table(phenotype)
pheno_color = c("darkolivegreen1", "darkolivegreen2", "gold3", "darkred")
#pheno_color = c("gold3", "darkred")
pheno_type = c("CN", "Obese", "NAFL", "NASH")

#
idx.c1 = grep("healthy", phenotype)
idx.c2 = grep("obese", phenotype)
idx.d1 = grep("NAFLD", phenotype)
idx.d2 = grep("NASH", phenotype)
phenotype[idx.c1] = pheno_type[1]
phenotype[idx.c2] = pheno_type[2]
phenotype[idx.d1] = pheno_type[3]
phenotype[idx.d2] = pheno_type[4]
names(pheno_color) = pheno_type

#
rowSds = rowSds(as.matrix(expr))
cut_off = quantile(rowSds, probs = 0.8)
expr = expr[rowSds > cut_off, ]

#
expr_raw = expr
expr = t(expr)
expr = scale(expr)
expr = t(expr)
quantile(expr)
expr[expr >= 5] = 5
expr[expr <= -5] = -5

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
pdf(figure_storage, height = 6, width = 10)

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
