library(data.table)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(scales)
library(reshape2)
library(plyr) 
library(grid)
library(gridExtra)
library(lattice)
library(patchwork)
library(devtools)
library(ggbiplot)
library(RColorBrewer)
library(vioplot)
library(ggridges)
library("stringr")

datadir="/Users/oushujun/Google Drive/study/Maize research/NAM/TE/bigNAM supplements/"
setwd(datadir)

# read data
fam = read.table('NAM.EDTA1.9.0.MTEC02052020.TE.v1.1.anno.sum.fam.bp', header=T)
class = read.table('NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.TEfam.list', header=F)
colnames(class) = c('TE', "Supfam") #add column header

# fix genome names
fam$TE_fam = sub("Oh7b", "Oh7B", fam$TE_fam)
names(fam)[names(fam) == 'Oh7b'] = 'Oh7B'
str(fam)

#remove AB10
drops <- c("B73Ab10","B73_AB10")
fam = fam[ , !(names(fam) %in% drops)]
colnames(fam) #check if Ab10 is removed

# transpose, and add back col and row names
tr_fam <- transpose(fam[-1]) # remove the first row, $TE_fam 
rownames(tr_fam) = colnames(fam[-1])
colnames(tr_fam) = as.matrix(transpose(fam[1]))[1,]
fam = tr_fam
fam = fam[,apply(fam, 2, mean) != 0] #remove all 0 rows
fam = fam/1000000 # convert bp to mb
fam = fam[-grep(":", colnames(fam),)] # remove rare TEs

## count family occurrance
levels=c("0")
non_0_count <- 26 - sapply(levels, function(x)colSums(fam=="0")) #count occurrences of 0 in each column

## add family and superfamily info and mean family size info
Fam_count = data.frame(count=non_0_count, TE=colnames(fam),
                       Supfam=class[match(colnames(fam), class$TE),]$Supfam,
                       Size=colSums(fam)/non_0_count)
colnames(Fam_count) = c('count', 'TE', 'Supfam', 'Size')


## filter out nonTE categories
Fam_count = Fam_count[which(Fam_count$Supfam!="centromeric_repeat" & Fam_count$Supfam!="knob" &
                              Fam_count$Supfam!="low_complexity" & Fam_count$Supfam!="rDNA_intergenic_spacer_element" & 
                              Fam_count$Supfam!="subtelomere"),]
Fam_count$Supfam = factor(Fam_count$Supfam) #remove unused factor levels
levels(Fam_count$Supfam)
str(Fam_count)

## replace LINEs to nonLTR
nonLTR = c("L1_LINE_retrotransposon", "RTE_LINE_retrotransposon", "LINE_element")
Fam_count$Supfam <- as.character(Fam_count$Supfam)
Fam_count$Supfam <- factor(with(Fam_count, replace(Supfam, Supfam %in% nonLTR, "nonLTR")))
levels(Fam_count$Supfam)

## Total number of TE families
length(Fam_count$count)

## Number of families that are present in all 26 NAM parents
Fam_count$TE = as.factor(Fam_count$TE)
sum = length(levels(Fam_count$TE[-grep(":", Fam_count$TE)])) #sum of TEs used in the NAM.EDTA.clean library
table(Fam_count$count)["26"]
table(Fam_count$count)["26"]/27228
table(Fam_count$count)["26"]/length(Fam_count$count)
table(Fam_count$count)["26"]/sum

## more than 20 times
sum(table(Fam_count$count)[20:26])
sum(table(Fam_count$count)[20:26])/27228
sum(table(Fam_count$count)[20:26])/length(Fam_count$count)
sum(table(Fam_count$count)[20:26])/sum

## Number of families that are only present in one NAM parent
table(Fam_count$count)["1"]
table(Fam_count$count)["1"]/27228
table(Fam_count$count)["1"]/length(Fam_count$count)
table(Fam_count$count)["1"]/sum
table(Fam_count[-grep(":", Fam_count$TE),]$count)["1"]/sum
table(Fam_count[-grep(":", Fam_count$TE),]$count)["1"]/27228

## other summaries
sum(table(Fam_count$count)[24:25])/27228
sum(table(Fam_count$count)[2:23])/27228
sum(table(Fam_count$count)[24:26])/27228

## histogram of freq
Fam_count$count = as.numeric(Fam_count$count)

Fam_count_freq_plot_all = ggplot(Fam_count, aes(count)) +
  geom_histogram(binwidth=1, center=0.5) +
  #  facet_wrap(~Supfam, scales = "free_y", ncol = 5) +
  labs(x ="Frequency (n = 26)", y = "Number of TE family") +
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.title.x = element_text(size=18, face="bold"),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_blank())
Fam_count_freq_plot_all

# save it to a pdf
pdf("Suppl. Fig. S5.pdf",width=6,height=10,pointsize=12, paper='special')
NAM.chr1_p
dev.off()