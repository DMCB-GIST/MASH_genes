rm(list = ls())
library(data.table)
setwd("c:/Dropbox/ddasic/GWAS/")

#
summary = fread("GCST90044186_buildGRCh37.tsv.gz")
summary = data.frame(summary)

#
title = "GCST90044186_LC"
colnames(summary)
snp_info = summary[, c("variant_id", "chromosome", "base_pair_location")]
colnames(snp_info) = c("snp", "chr", "pos")

#
data = data.frame(
  chr = snp_info$chr,
  pos = snp_info$pos
)
data$pos = as.numeric(data$pos)
data$chr = as.numeric(data$chr)

# cM 단위로 윈도우 크기 설정 (예: 1cM)
cM_window <- 1  # cM 단위로 변경 가능

# 각 크로모좀에 대해 cM 윈도우로 데이터 분할 및 SNP 밀도 계산
data <- data %>%
  mutate(cM_position = pos / 1e7,  # 예제에서는 pos를 cM으로 변환 (실제 cM 변환 필요)
         window = cut(cM_position, breaks = seq(0, (max(cM_position)+1), by = cM_window), labels = FALSE))%>%
  group_by(chr, window) %>%
  summarise(SNP_count = n(), .groups = 'drop')

#
library(ggplot2)
hist(data$SNP_count)
breaks = c(10, 10000, 20000, 30000, 40000, 50000, 60000)
label = c(breaks, paste0(">", max(breaks)))
breaks = c(-Inf, breaks, Inf)
data <- data %>%
  mutate(SNP_bin = cut(SNP_count, 
                       breaks = breaks, 
                       labels = label))
#colors = colorRampPalette(c("#fdf498", "#7bc043", "#ff914d", "red"))(length(breaks)-1)
colors = colorRampPalette(c("#f7f7f7", "#dfe3ee", "#8b9dc3", "#3b5998"))(length(breaks)-1)

#
# data = rbind(c(1, 1, 1, breaks[2]), data)
# data = rbind(c(1, 2, 1, breaks[3]), data)
# data = rbind(c(1, 3, 1, breaks[4]), data)
# data = rbind(c(1, 4, 1, breaks[5]), data)
data = data[!is.na(data$chr), ]

p = ggplot(data, aes(x = window, y = factor(chr, levels = 22:1), fill = SNP_bin)) +
  geom_tile(color = "white") +
  #scale_fill_gradientn(colours = , na.value = "white") +
  #scale_fill_gradientn(colors = colors, breaks = breaks, na.value = "white") +
  scale_fill_manual(values = colors, na.value = "white")+
  #scale_x_continuous(expand = c(0, 0), limits = c(1500, 1750)) +
  labs(x = "Genetic Position (cM)", y = "Chromosome", fill = "SNP Count",
       title = title) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing = unit(0, "lines"),
        panel.grid.minor.x = element_blank())
wd = "c:/Dropbox/ddasic/05.fri/01.NASH_gene/data/"
storage = paste0(wd, "FS.", title, "_density.png")
ggsave(storage, width = 6, height = 4)

#
summary$variant_id[1:100]
# source("c:/Dropbox/ddasic/statistics/meta_analysis_making_object.r")
# files = list.files(pattern = "deseq2")
# if(length(files) > 0){
#   for(i in 1:length(files)){
#     summary = read.table(files[i])
#     summary$marker = rownames(summary)
#     summary$effect = summary$log2FoldChange
#     summary$stderr = summary$lfcSE
#     summary = summary[, c("marker", "effect", "stderr")]
#     summary = summary[complete.cases(summary), ]
#     file_storage = gsub("deseq2", "metal", files[i])
#     write.table(summary, file_storage, sep = "\t", row.names = F, quote = F)
#   }
summary$marker = summary$variant_id
summary$chr = summary$chromosome  
summary$pos = summary$base_pair_location
summary$effect = summary$beta
summary$stderr = summary$standard_error
summary$p = summary$p_value
summary = summary[, c("marker", "effect", "stderr", "p", "chr", "pos")]
write.table(summary, "GCST90044186_LC_processed.txt", sep = "\t", row.names = F, quote = F)
