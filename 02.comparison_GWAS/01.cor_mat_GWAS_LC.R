rm(list = ls())
wd = "c:/Dropbox/ddasic/05.fri/01.NASH_gene/data/"
setwd(wd)
library(data.table)

#
files = c("GCST90018606_LC_metal.txt", 
          "GCST90018826_LC_metal.txt",
          "GCST90044186_LC_metal.txt",
          "GCST90319877_LC_metal.txt",
          "GCST90436370_LC_metal.txt")

#
names = c("GCST90018606", 
              "GCST90018826",
              "GCST90044186",
              "GCST90319877",
              "GCST90436370")

#
df = NULL
for(i in 1:length(files)){
  res.i = fread(files[i])
  res.i = data.frame(res.i)
  t.i = res.i$effect/res.i$stderr
  for(j in 1:length(files)){
    if(i < j){
      res.j = fread(files[j])
      idx1 = match(res.i$marker, res.j$marker)
      res.j = res.j[idx1, ]
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
ggsave("FS.GWAS_comparison_LC.pdf", width = 6, height = 6)
