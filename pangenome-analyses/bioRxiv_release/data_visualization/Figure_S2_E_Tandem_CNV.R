library(plotrix)
library(tidyr)
library(matrixStats)
library(tidyr)

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


visualization_tandem_CNV <- read.csv("visualization_tandem_CNV_copy.csv")


visualization_tandem_CNV$bin <-ifelse(visualization_tandem_CNV$mean>=2 & visualization_tandem_CNV$mean <2.5,rr2<-"2",
                             ifelse(visualization_tandem_CNV$mean>=2.5 & visualization_tandem_CNV$mean<3.5,rr2<-"3",
                                    ifelse(visualization_tandem_CNV$mean>=3.5 & visualization_tandem_CNV$mean <4.5,rr2<-"4",
                                           ifelse(visualization_tandem_CNV$mean>=4.5 & visualization_tandem_CNV$mean <5.5,rr2<-"5",
                                                  ifelse(visualization_tandem_CNV$mean>=5.5 & visualization_tandem_CNV$mean <6.5,rr2<-"6",
                                                         ifelse(visualization_tandem_CNV$mean>=6.5 & visualization_tandem_CNV$mean<7.5,rr2<-"7",
                                                                ifelse(visualization_tandem_CNV$mean>=7.5 & visualization_tandem_CNV$mean <8.5,rr2<-"8",
                                                                       ifelse(visualization_tandem_CNV$mean>=8.5 & visualization_tandem_CNV$mean <9.5,rr2<-"9",
                                                                                     rr2<-">=10"))))))))
visulization_tandem_mean <- as.data.frame(table(visualization_tandem_CNV$bin))

visulization_tandem_mean$Var1 <- factor(visulization_tandem_mean$Var1, levels=c("2", "3","4","5","6","7","8","9",">=10"))


ggplot(visulization_tandem_mean, aes(x = Var1, y = Freq)) + theme_bw() + geom_bar(stat = "identity") + geom_text(aes(label=Freq),vjust=-0.5) +
  ylab("Number of Putative Tandem Duplicates") + xlab("Putative Tandem Duplicates Mean Copy Number") + theme_classic() + 
  theme(text = element_text(size = 12),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))



