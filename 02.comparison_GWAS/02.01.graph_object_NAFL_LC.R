rm(list = ls())
wd = "c:/Dropbox/ddasic/05.fri/01.NASH_gene/data/"
setwd(wd)
library(data.table)

#

files = c("GCST90011885_NAFL_metal.txt", 
          "GCST90054782_NAFL_metal.txt", 
          "GCST90091033_NAFL_metal.txt",
          "GCST90275041_NAFL_metal.txt",
          "GCST90018606_LC_metal.txt", 
          "GCST90018826_LC_metal.txt",
          "GCST90044186_LC_metal.txt",
          "GCST90319877_LC_metal.txt",
          "GCST90436370_LC_metal.txt")

#
mat = matrix(0, nrow = length(files), ncol = length(files))
colnames(mat) = rownames(mat) = files
for(i in 1:length(files)){
  res.i = fread(files[i])
  res.i = data.frame(res.i)
  for(j in 1:length(files)){
    if(i < j){
      res.j = fread(files[j])
      common_snp = intersect(res.i$marker, res.j$marker)
      idx1 = match(common_snp, res.i$marker)
      res.i1 = res.i[idx1, ]
      idx1 = match(common_snp, res.j$marker)
      res.j = res.j[idx1, ]
      t.i = res.i1$effect/res.i1$stderr
      t.j = res.j$effect/res.j$stderr
      t = data.frame(i = t.i, j = t.j)
      t = t[complete.cases(t), ]
      cor.test = cor.test(t[, 1], t[, 2])
      pcc = cor.test$estimate
      mat[i, j] = pcc
      mat[j, i] = pcc
    }
  }
}

save(file = "GWAS_network_obj.rdata", mat)
