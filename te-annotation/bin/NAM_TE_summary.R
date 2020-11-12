# Shujun Ou (shujun.ou.1@gmail.com) 
# Sept. 2020

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

datadir="/Users/oushujun/Google Drive/study/Maize research/NAM/TE"
setwd(datadir)


########################
#### basic settings ####
########################

NAM_ID = c("B73", "B97", "Ky21", "M162W", "Ms71", "Oh7b", "Oh43", "M37W", "Mo18W", "Tx303",
           "HP301", "Il14H", "P39", "CML52", "CML69", "CML103", "CML228", "CML247", "CML277", 
           "CML322", "CML333", "Ki3", "Ki11", "NC350", "NC358", "Tzi8")
NAM_colors = c("B73"="goldenrod1", "B97"="royalblue", "CML103"="limegreen", 
               "CML228"="limegreen", "CML247"="limegreen", "CML277"="limegreen", 
               "CML322"="limegreen", "CML333"="limegreen", "CML52"="limegreen", 
               "CML69"="limegreen", "HP301"="orchid", "Il14H"="orangered", 
               "Ki11"="limegreen", "Ki3"="limegreen", "Ky21"="royalblue", 
               "M162W"="royalblue", "M37W"="gray47", "Ms71"="royalblue", 
               "Mo18W"="gray47", "NC350"="limegreen", "NC358"="limegreen", 
               "Oh43"="royalblue", "Oh7b"="royalblue", "P39"="orangered", 
               "Tx303"="gray47", "Tzi8"="limegreen")
#reorder NAM_colors based on NAM_ID
NAM_colors = NAM_colors[order(match(names(NAM_colors), NAM_ID))]
TE_colors = c('gray', rev(brewer.pal(11,'RdYlBu')[1:5]), 
              brewer.pal(11,'RdYlBu')[6], 
              rev(brewer.pal(11,'RdBu')[7:10]))

show_col(TE_colors)

#############
#### LAI ####
#############
LAI = read.table('NAM_LAI_v2.8.txt', header=T)

# remove the B73_ab10 column
LAI = LAI[-1,]
summary(LAI$LAI[4:29])
sd(LAI$LAI[4:29])
str(LAI)

# change the order of Genomes
LAI_ID = c("W22v2", "Mo17v1", "B73v4", "B73v5", NAM_ID[-1])
LAI_colors = c("W22v2"='gray', "Mo17v1"='gray', "B73v4"='gray', "B73v5"='goldenrod1', NAM_colors[-1])
LAI_colors = LAI_colors[order(match(names(LAI_colors), LAI_ID))]
LAI$Genome = factor(LAI$Genome, levels=LAI_ID)

#plot
LAI_plot = LAI %>%
  ggplot(aes(x = Genome, y = LAI, fill = Genome)) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values=LAI_colors) +
  labs(x =" ", y = "LTR Assembly Index (LAI)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=16), 
        axis.text.x = element_text(face="bold", size=16, angle=35, vjust=1, hjust=0.95)) +
  theme(legend.position="")
LAI_plot

#output
pdf("NAM.TE.LAI.pdf", width=11,height=6,pointsize=12, paper='special')
LAI_plot
dev.off()


############################
#### TE size stack plot ####
############################
## total TE bp in assembly
repeats = read.table('NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.bp.txt', header=T)

# add more categories
repeats$DNA = repeats$helitron + repeats$CACTA + repeats$Mutator + repeats$PIF_Harbinger +
  repeats$Tc1_Mariner + repeats$hAT
repeats$nonLTR = repeats$L1_LINE + repeats$LINE_element + repeats$RTE_LINE
repeats$LTR = repeats$Copia + repeats$Gypsy + repeats$LTR_unknown
repeats$nonTE = repeats$centromeric_repeat + repeats$knob + repeats$low_complexity + 
  repeats$rDNA_spacer + repeats$subtelomere
repeats$TE = repeats$nonLTR + repeats$LTR + repeats$DNA

#format genome name
repeats$Genome = sub("_bp", "", repeats$Genome, fixed=TRUE)
head(repeats)
str(repeats)

## total intact TE bp in assembly
intact = read.table('NAM.EDTA1.9.0.MTEC02052020.intact.sum.txt', header=T)

#remove Ab10 columns
repeats = repeats[-c(1,2),]
intact = intact[-c(1,2),]
intact$intact = intact$LTR + intact$TIR + intact$Helitron
intact$fragmented = repeats$TE - intact$intact
str(intact)
intact$Genome

#read in genome size data
genome_size = read.table('Genome_size.txt', header=T)
genome_size = genome_size[order(match(genome_size$Genome, repeats$Genome)),]

#summarize intact and fragmented percentage
repeats_pcnt = data.frame(repeats[,-1]/10000000/genome_size$Size, Genome=repeats$Genome)
repeats_bp = data.frame(repeats[,-1]/1000000, Genome=repeats$Genome)
intact_pcnt = data.frame(intact[,-1]/10000000/genome_size$Size, Genome=intact$Genome)
intact_bp = data.frame(intact[,-1]/1000000, Genome=intact$Genome)
intact_pcntTE = data.frame(intact[,-1]*100/(intact$intact+intact$fragmented), Genome=intact$Genome) #intact/TE, frag/TE
summary(intact_pcntTE$intact)
summary(intact_pcntTE$fragmented)

# change the order of Genomes
repeats_pcnt$Genome = factor(repeats_pcnt$Genome, levels=NAM_ID)
repeats_bp$Genome = factor(repeats_bp$Genome, levels=NAM_ID)


# plot the %intact and %homo in %TE
intact_pcntTE_final = intact_pcntTE %>% gather(variable, value, intact, fragmented)
intact_pcntTE_final$variable = as.factor(intact_pcntTE_final$variable) #convert to character to factor
intact_pcntTE_final$variable = factor(intact_pcntTE_final$variable, 
                                      levels = c("intact", "fragmented"))

# plot the %intact and %frag to the total TE size
intact_pcntTE_plot = intact_pcntTE_final %>%
  ggplot(aes(x = Genome, y = value, fill = variable)) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values=c("#92C5DE", "gray")) +
  labs(x =" ", y = "TE contents (%)",size=14, face="bold") +
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(color = NAM_colors, face="bold", size=16, angle=35, vjust=1, hjust=0.95),
        legend.text = element_text(size=14),
        legend.title=element_blank()) +
  theme(legend.position="top")
intact_pcntTE_plot


# gather plotting variables and change the variable plotting order
repeats_final = repeats_bp %>% 
  gather(variable, value, hAT, CACTA, PIF_Harbinger, Mutator, Tc1_Mariner,  
         helitron, nonLTR, Copia, Gypsy, LTR_unknown) #remove nonTE
repeats_final$variable = as.factor(repeats_final$variable) #convert to character to factor
repeats_final$variable = factor(repeats_final$variable, 
                                levels = c("hAT", "CACTA", "PIF_Harbinger", "Mutator", 
                                           "Tc1_Mariner", "helitron", "nonLTR", 
                                           "LTR_unknown", "Gypsy", "Copia"))

# plot the total Mb of TEs
repeats_Mb_plot = repeats_final %>%
  ggplot(aes(x = Genome, y = value, fill = variable)) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values=TE_colors[-1]) + #remove nonTE
  labs(x =" ", y = "Repeat contents (Mb)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(color = NAM_colors, face="bold", size=16, angle=35, vjust=1, hjust=0.95),
        legend.text = element_text(size=14),
        legend.title=element_blank()) +
  theme(legend.position="top")
repeats_Mb_plot

# add percent to B73 as the second axis
B73_size = 2178268108/1000000
repeats_Mb_plot_B73pcnt = 
  repeats_Mb_plot + 
  scale_y_continuous(sec.axis = sec_axis(~ . *100/B73_size, name="Size to B73 (%)")) # +
theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 50)))
repeats_Mb_plot_B73pcnt

## test if tropical lines has more TEs than non-tropical lines
str(repeats_bp)

tropical = c("CML103", "CML228", "CML247", "CML277", "CML322", "CML333", "CML52", "CML69", "Ki3", "Ki11", "NC350", "NC358", "Tzi8")
non_tropical = c("B73", "B97", "Ky21", "M162W", "MS71", "Oh43", "Oh7b", "M37W", "Mo18W", "Tx303", "HP301", "P39", "IL14H")
repeat_trop = repeats_bp[repeats_bp$Genome %in% tropical,]
repeat_notrop = repeats_bp[repeats_bp$Genome %in% non_tropical,]

#For loop
testresults <- vector("list", ncol(repeats_bp)-1)
for (j in seq(ncol(repeats_bp)-1)){
  testresults[[j]] = t.test(repeat_trop[,j], repeat_notrop[,j])
}
testresults
colnames(repeat_trop)

#These are significant: 
#trop significantly larger than non-trop
#gypsy (2, p=0.000739), CACTA (4, p=0.01612), Mutator (5, p=0.01669), PIF (6, p=0.04081),
#Tc1 (7, p=0.03369), Total (18, p=0.0003768), LTR (21, p=0.0002125), TE (23, p=0.0004438)
#real difference: trop >> non-trop
#gypsy (2), CACTA (4), Mutator (5), PIF (6)
#real difference: trop << non-trop
#Tc1 (7)

#output
pdf("NAM.TE.summary2.pdf", width=11,height=6,pointsize=12, paper='special')
intact_pcntTE_plot
repeats_Mb_plot_B73pcnt
dev.off()


######################################
#### Read and process family data ####
######################################
fam = read.table('NAM.EDTA1.9.0.MTEC02052020.TE.v1.1.anno.sum.fam.bp', header=T)
class = read.table('NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.TEfam.list', header=F)
colnames(class) = c('TE', "Supfam") #add column header

#remove AB10
drops <- c("B73Ab10","B73_AB10")
fam = fam[ , !(names(fam) %in% drops)]
colnames(fam) #check if Ab10 is removed

# transpose, and add back col and row names
library(data.table)
tr_fam <- transpose(fam[-1]) # remove the first row, $TE_fam 
rownames(tr_fam) = colnames(fam[-1])
colnames(tr_fam) = as.matrix(transpose(fam[1]))[1,]

#str(tr_fam)
fam = tr_fam
fam = fam[,apply(fam, 2, mean) != 0] #remove all 0 rows

# convert bp to mb
fam = fam/1000000

# convert bp to TE space percentage
fam_pcnt = fam/rowSums(fam)


############################
#### Pan TE family stat ####
############################
## Remove entries that are not in the pan-TE library
fam = fam[-grep(":", colnames(fam),)]

## count family occurrance
str(fam)
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

## histogram of freq
plotorder <- c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", "LTR_retrotransposon", 
               "nonLTR", "helitron", "CACTA_TIR_transposon", "hAT_TIR_transposon", "Mutator_TIR_transposon", 
               "PIF_Harbinger_TIR_transposon", "Tc1_Mariner_TIR_transposon")
Fam_count <- arrange(transform(Fam_count,
                               Supfam=factor(Supfam,levels=plotorder)),Supfam)

Fam_count$count = as.numeric(Fam_count$count)
Fam_count_freq_plot = ggplot(Fam_count, aes(count)) +
  geom_histogram(binwidth=1) +
  facet_wrap(~Supfam, scales = "free_y", ncol = 5) +
  labs(x ="Frequency (n = 26)", y = "Number of TE family") +
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.title.x = element_text(size=18, face="bold"),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_blank())
Fam_count_freq_plot

Fam_count_freq_plot_all = ggplot(Fam_count, aes(count)) +
  geom_histogram(binwidth=1) +
  #  facet_wrap(~Supfam, scales = "free_y", ncol = 5) +
  labs(x ="Frequency (n = 26)", y = "Number of TE family") +
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.title.x = element_text(size=18, face="bold"),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_blank())
Fam_count_freq_plot_all

#output
pdf("NAM.TE.family2.pdf", width=11,height=6,pointsize=12, paper='special')
Fam_count_freq_plot
Fam_count_freq_plot_all
dev.off()
