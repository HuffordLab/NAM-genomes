tandem_CNV <- read.csv('~/Desktop/NAM_PAN_GENOME/pan_genome_nov 2/QC_set/tandem_QC/tanden_gene_type_CNV_mean.csv')

# replace 0 and 1 copy into missing data 
library(naniar)
library(tidyr)
tandem_CNV[2:28] %>% replace_with_na_all(condition = ~.x == 1) %>% replace_with_na_all(condition = ~.x == 0) %>% write.csv(file = '~/Desktop/tanden_cnv_mean.csv')



# modify the query gene by copying from the orignial matrix 
library(ggplot2)
library(plotrix)
library(tidyr)
library(matrixStats)
visualization_tandem_CNV <- read.csv('~/Desktop/NAM_PAN_GENOME/pan_genome_nov 2/QC_set/tandem_QC/tanden_cnv_with_mean_value.csv')
# mean of the mean
mean(visualization_tandem_CNV$mean)
#[1] 2.199297
# sd of the mean
std.error(visualization_tandem_CNV$mean)
# [1] 0.007200943


# whole set 
table(visualization_tandem_CNV$mean)
whole_set <- ggplot(visualization_tandem_CNV,aes(mean)) + stat_bin(binwidth=1) + #stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5) + 
  scale_x_continuous(breaks = seq(0, 60, by = 5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),text = element_text(size = 20)) + 
  #scale_y_continuous(breaks = seq(0, 10, by = 1))  +
  ylab("Number of Putative Tandem Duplicates") + xlab("Putative Tandem Duplicates Copy Number Mean (2-60)") + theme_classic()


#subset for zoom in 
copy_variation_10_and_below <- subset(visualization_tandem_CNV,visualization_tandem_CNV$mean <=10)
table(copy_variation_10_and_below$mean)
# creating grobs to add in 
bp<- ggplotGrob(ggplot(copy_variation_10_and_below,aes(mean)) + stat_bin(binwidth=0.5) + #stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5,angle =0,colour="red",size=3) + 
                  scale_x_continuous(breaks = seq(0, 10, by = 1)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),text = element_text(size = 12)) + theme_classic() +
                  ylab("Number of Putative Tandem Duplicates") +  xlab("Putative Tandem Duplicates Copy Number Mean (2-10)"))


tanden_cnv <- whole_set + 
  annotation_custom(
    grob = bp,
    xmin = 10,
    xmax = 55,
    ymin = 1000,
    ymax = 11000
  ) 
tanden_cnv + theme(text = element_text(size = 14))
