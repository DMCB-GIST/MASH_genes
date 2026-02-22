rm(list = ls())
wd = "c:/Dropbox/ddasic/05.fri/01.NASH_gene/data/"
setwd(wd)

#
files = list.files(pattern = "processed")

#
names = gsub("_processed.txt", "", files)

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
df = NULL
for(i in 1:length(files)){
  res.i = fread(files[i])
  res.i = data.frame(res.i)
  res.i = res.i[res.i$p < 0.2, ]
  t.i = res.i$effect/res.i$stderr
  for(j in 1:length(files)){
    if(i < j){
      res.j = fread(files[j])
      res.j = res.j[res.j$p < 0.2, ]
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
df_02 = df

#
df = NULL
p = 0.15
for(i in 1:length(files)){
  res.i = fread(files[i])
  res.i = data.frame(res.i)
  res.i = res.i[res.i$p < p, ]
  t.i = res.i$effect/res.i$stderr
  for(j in 1:length(files)){
    if(i < j){
      res.j = fread(files[j])
      res.j = res.j[res.j$p < p, ]
      idx1 = match(res.i$marker, res.j$marker)
      res.j = res.j[idx1, ]
      t.j = res.j$effect/res.j$stderr
      t = data.frame(i = t.i, j = t.j)
      t = t[complete.cases(t), ]
      cor.test = try(cor.test(t[, 1], t[, 2]))
      if(length(cor.test) != 1){
        pcc = cor.test$estimate
        p = cor.test$p.value
        df.i = data.frame(group = paste0(names[i], "_", names[j]), pcc = pcc, p = p, n = nrow(t))
        df = rbind(df, df.i)
      } else {
        df.i = data.frame(group = paste0(names[i], "_", names[j]), pcc = NA, p = NA, n = NA)
      }
    }
  }
}
df_015 = df

#
df_all$fill = "All"
df_all = df_all[order(-1*df$pcc), ]
df_02$fill = "p02"

#
df = rbind(df_all, df_02)
df$fill = as.factor(df$fill)

#
library(ggplot2)
df$group = factor(df$group, levels = df_all$group)
ggplot(df, aes(x = group, y = pcc, fill = fill))+
  geom_bar(stat = "identity", position = position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
ggsave("FS.GWAS_comparison.pdf", width = 6, height = 6)
