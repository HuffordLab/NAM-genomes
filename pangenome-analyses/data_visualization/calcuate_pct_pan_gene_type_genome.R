plot_gene_type <- read.csv("summary_gene_type_fmt_final.csv")

library(plyr)
summ <- ddply(plot_gene_type, .(NAM, Gene_Type), summarize, Sum_gene=sum(Gene_Count))
gene_type_percentage <- ddply(summ, .(NAM), mutate, gene_pct = Sum_gene / sum(Sum_gene) * 100)

write.csv(as.data.frame(gene_type_percentage),file="gene_type_percentage.csv")

# reshape the matrix only using percentage 
# calculate percentage and standard error 
library(plotrix)
se_calculation <- read.csv("standard_error_final.csv")
se_calculation
se_calculation

core_mean <- mean(se_calculation$core)
core_mean
std.error(se_calculation$core)

softcore_mean <-mean(se_calculation$softcore)
softcore_mean
std.error(se_calculation$softcore)

dis_mean <- mean(se_calculation$dispensable)
dis_mean
std.error(se_calculation$dispensable)

pri_mean <- mean(se_calculation$private)
pri_mean
std.error(se_calculation$private)
core_mean + softcore_mean+dis_mean + pri_mean


