rm(list = ls())
wd = "c:/Dropbox/ddasic/05.fri/01.NASH_gene/data/"
setwd(wd)
library(data.table)

#
files = c("GCST90011885_NAFL_processed.txt", 
          "GCST90054782_NAFL_processed.txt", 
          "GCST90091033_NAFL_processed.txt",
          "GCST90271622_NAFL_processed.txt",
          "GCST90275041_NAFL_processed.txt")

#
names = c("GCST90011885_NAFL", 
              "GCST90054782_NAFL",
              "GCST90091033_NAFL",
              "GCST90271622_NAFL",
              "GCST90275041_NAFL")

#
df = NULL
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
      p = cor.test$p.value
      df.i = data.frame(group = paste0(names[i], "_", names[j]), pcc = pcc, p = p, n = nrow(t))
      df = rbind(df, df.i)
    }
  }
}
df_all = df

#
library(ggplot2)
df$group = factor(df$group, levels = df_all$group)
ggplot(df, aes(x = group, y = pcc))+
  geom_bar(stat = "identity", position = position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
  ylab("PCC")
ggsave("FS.GWAS_comparison_NAFL.pdf", width = 6, height = 6)
