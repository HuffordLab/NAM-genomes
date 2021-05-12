## Summary of TEs in NAM genomes
## Shujun Ou (shujun.ou.1@gmail.com)
## 01/05/2021

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


########################
#### basic settings ####
########################

NAM_ID = c("B73", "B97", "Ky21", "M162W", "Ms71", "Oh43", "Oh7B", "M37W", "Mo18W", "Tx303",
           "HP301", "P39", "Il14H", "CML52", "CML69", "CML103", "CML228", "CML247", "CML277", 
           "CML322", "CML333", "Ki3", "Ki11", "NC350", "NC358", "Tzi8")

NAM_colors = c("B73"="goldenrod1", "B97"="royalblue", "CML103"="limegreen", 
               "CML228"="limegreen", "CML247"="limegreen", "CML277"="limegreen", 
               "CML322"="limegreen", "CML333"="limegreen", "CML52"="limegreen", 
               "CML69"="limegreen", "HP301"="orchid", "Il14H"="orangered", 
               "Ki11"="limegreen", "Ki3"="limegreen", "Ky21"="royalblue", 
               "M162W"="royalblue", "M37W"="gray47", "Ms71"="royalblue", 
               "Mo18W"="gray47", "NC350"="limegreen", "NC358"="limegreen", 
               "Oh43"="royalblue", "Oh7B"="royalblue", "P39"="orangered", 
               "Tx303"="gray47", "Tzi8"="limegreen")
#reorder NAM_colors based on NAM_ID
NAM_colors = NAM_colors[order(match(names(NAM_colors), NAM_ID))]
TE_colors = c('gray', rev(brewer.pal(11,'RdYlBu')[1:5]), 
              brewer.pal(11,'RdYlBu')[6], 
              rev(brewer.pal(11,'RdBu')[7:10]))

show_col(TE_colors)

############################
#### Suppl. Fig. S4 ####
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

#remove Ab10 columns
repeats = repeats[-c(1,2),]

# fix genome names
repeats$Genome = sub("Oh7b", "Oh7B", repeats$Genome)

# change the order of Genomes
repeats_bp = data.frame(repeats[,-1]/1000000, Genome=repeats$Genome)
repeats_bp$Genome = factor(repeats_bp$Genome, levels=NAM_ID)

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
  scale_y_continuous(sec.axis = sec_axis(~ . *100/B73_size, name="Size to B73 (%)"))
repeats_Mb_plot_B73pcnt

#output
pdf("Suppl. Fig. S4.pdf", width=11,height=6,pointsize=12, paper='special')
  repeats_Mb_plot_B73pcnt
dev.off()


