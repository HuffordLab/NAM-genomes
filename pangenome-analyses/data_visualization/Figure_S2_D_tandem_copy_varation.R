library(naniar)
library(tidyr)
library(ggplot2)
library(plotrix)
library(tidyr)
library(matrixStats)

# Tandem duplicate copy number mean 

number_matrix <- read.csv("numeric_pan_matrix.csv",header= TRUE) 
tandem_gene_id <- read.csv(file = "tandem_duplicate_ID.txt",sep="\t")

tandem_copy_matrix <- left_join(tandem_gene_id,number_matrix[,-1])


# replace 0 and 1 copy into missing data 
visualization_tandem_CNV <- tandem_copy_matrix[2:28] %>% replace_with_na_all(condition = ~.x == 1) %>% replace_with_na_all(condition = ~.x == 0) 
visualization_tandem_CNV$mean<- rowMeans(tandem_convert[1:26],na.rm = TRUE)
# count how many NA per row, if NA ==25, meaning this tandem is private to one genotype 

visualization_tandem_CNV$NA_count 
visualization_tandem_CNV$na_count <- apply(visualization_tandem_CNV[1:26], 1, function(x) sum(is.na(x)))
table(visualization_tandem_CNV$na_count)

mean(visualization_tandem_CNV$mean)
#[1] 2.196584 
# sd of the mean
std.error(visualization_tandem_CNV$mean)
# [1] 0.007016222

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



# response to reviewer 

cnv <- as.matrix(visualization_tandem_CNV[,1:26])
cnv[is.na(cnv)] <- 0

copy_number_range <- as.data.frame(rowRanges(cnv, na.rm = TRUE, dims = 1, n = NULL))
colnames(copy_number_range) <- c("Min_copy","Max_copy")

same_copy = copy_number_range %>% filter(Min_copy == Max_copy)
# 237 same copy

