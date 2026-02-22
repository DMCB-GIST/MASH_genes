rm(list = ls())
wd = "c:/Dropbox/ddasic/05.fri/01.NASH_gene/data/"
setwd(wd)
library(data.table)

#
load("GWAS_network_obj.rdata")
#source("c:/Dropbox/ddasic/unsupervied_analysis/network_analysis.r")
#source("c:/Dropbox/ddasic/unsupervied_analysis/network_plot.r")

#
library(igraph)
library(ggplot2)

#
col_n = rownames(mat)
col_n = gsub("_metal.txt", "", col_n)
rownames(mat) = colnames(mat) = col_n
cor_mat = mat
# se = sqrt((1 - cor_mat^2)/(nrow(cor_mat)-2))
# cor_mat = cor_mat/se
# cor_mat[cor_mat >= 10] = 10
# cor_mat[cor_mat <= -10] = -10
# #hist(cor_mat)
# cor_mat = cor_mat/10

#
hist(cor_mat)

#
diag(cor_mat) = 0
cor_mat[is.na(cor_mat)] = 0

graph = graph.adjacency(cor_mat, weighted=TRUE, mode="lower")

#
cut_off_pcc = 0.03
g = delete.edges(graph, E(graph)[abs(weight) < cut_off_pcc])
g = simplify(g)
gsize(graph)
#node.shape#square#circle
V(g)$shape = "circle"

#degree
degree = degree(g, mode="all")
#hist(degree)
cut_off = c(1, 2)
idx1 = which(degree >= cut_off[1])
idx2 = which(degree >= cut_off[2])

#nodel.color
#V(g)$color = adjustcolor("#FADA2A", alpha.f = 1)
idx_n1 = grep("NAFL", names(V(g)))
idx_n2 = grep("LC", names(V(g)))
V(g)[idx_n1]$color = adjustcolor("#FADA2A", alpha.f = 1)
V(g)[idx_n2]$color = adjustcolor("brown3", alpha.f = 1)

#node.size
V(g)$size = 5
V(g)[degree > cut_off[1]]$size = 7
V(g)[degree > cut_off[2]]$size = 9

#lable.size
V(g)$label.cex = 1
#V(g)[degree > cut_off[1]]$label.cex = 0.
#V(g)[degree > cut_off[2]]$label.cex = 2

#label.color
V(g)$label.color = "grey30"

#
#E(g)$width = 2

#edge_color
idx1 = which(E(g)$weight < 0)
idx2 = which(E(g)$weight > 0)
E(g)[idx1]$color = adjustcolor("#5B89B4", alpha.f = 0.7)
E(g)[idx2]$color = adjustcolor("#CA7066", alpha.f = 0.7)

pdf("F.GWAS_network.pdf", width = 8, height = 8)
plot(g,
     #
     edge.width=abs(E(g)$weight)*5,
     edge.arrow.size = 0,
     #edge.curved = 0.5,
     #
     #layout=layout.fruchterman.reingold,
     #
     vertex.frame.color = NA,
     vertex.label.degree = 0.5,
     vertex.label.font=2, # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.shape="circle", # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
     vertex.label.dist = 1.2
)
dev.off()
