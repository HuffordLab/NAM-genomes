# check if ggplot2 is installed, if so, load it, 
# if not, install and load it
if("ggplot2" %in% rownames(installed.packages())){
  library(ggplot2)
} else {
  install.packages("ggplot2")
  library(ggplot2)
}
library('dplyr')
library(tidyr)
library(gridExtra)
library(grid)
library(lattice)
library(cowplot)
library(forcats)

datadir="/Users/oushujun/Google Drive/study/Maize research/NAM/TE"
setwd(datadir)

## plot NAM stats on Chr1
NAM.intact = read.table('NAM.EDTA.TEanno.structural.list', header=F)
NAM.frag = read.table('NAM.EDTA.TEanno.homology.list', header=F)
NAM.rmap = read.table('/Users/oushujun/Google Drive/study/Maize research/NAM/NAM-recombination/VCF_files_v4/NAM.crossover.chr1.pos', header=F)
NAM.intact.bp = read.table('NAM.EDTA.TEpass.sum.chr1.list', header=F)
NAM.all.bp = read.table('NAM.EDTA.TEall.sum.chr1.list', header=F)
colnames(NAM.frag) = c('chr', 'pos')
colnames(NAM.intact) = c('chr', 'pos')
colnames(NAM.rmap) = c('panel', 'chr', 'pos', 'marker')
colnames(NAM.intact.bp) = c("Chr","From","To","Copia","Gypsy","unknown","Total_LTR")
colnames(NAM.all.bp) = c("Chr","From","To","Copia","Gypsy","unknown","Total_LTR")
NAM.rmap$chr=paste("B73_chr", NAM.rmap$chr, sep="")
str(NAM.frag)
str(NAM.rmap)
head(NAM.rmap, 10)
str(NAM.intact.bp)

#estimate raw LAI
NAM.all.bp$intact.pcnt = NAM.intact.bp$Total_LTR/NAM.all.bp$Total_LTR

# plot chr1 of NAM 
NAM_intact_p <- ggplot() + 
  geom_density(data=NAM.intact, aes(x=pos/1e6, group=chr, color=chr), adjust=0.05, lwd=0.05, trim=T) + 
  ylab("Intact LTR Density") +
  theme(axis.title.y = element_text(size=9), legend.position = "none", axis.title.x=element_blank())
NAM_frag_p <- ggplot() + 
  geom_density(data=NAM.frag, aes(x=pos/1e6, group=chr, color=chr), adjust=0.05, lwd=0.05, trim=T) + 
  ylab("Fragmented LTR Density") +
  theme(axis.title.y = element_text(size=9), legend.position = "none", axis.title.x=element_blank())
NAM.chr1_total_p <- ggplot() + 
  geom_line(data=NAM.all.bp, aes(x = From / 1e6, y = Total_LTR * 100, color=Chr), lwd=0.05, position = "dodge") +
  ylab("Total LTR (%)") +
  theme(axis.title.y = element_text(size=9), legend.position = "none", axis.title.x=element_blank())
NAM.chr1_gypsy_p <- ggplot() + 
  geom_line(data=NAM.all.bp, aes(x = From / 1e6, y = Gypsy * 100, color=Chr), lwd=0.001, position = "dodge") +
  ylab("Total Gypsy LTR (%)") +
  theme(axis.title.y = element_text(size=9), legend.position = "none", axis.title.x=element_blank())
NAM.chr1_copia_p <- ggplot() + 
  geom_line(data=NAM.all.bp, aes(x = From / 1e6, y = Copia * 100, color=Chr), lwd=0.001, position = "dodge") +
  ylab("Total Copia LTR (%)") +
  theme(axis.title.y = element_text(size=9), legend.position = "none", axis.title.x=element_blank())
NAM.chr1_intact_adj_p <- ggplot() + 
  geom_line(data=NAM.all.bp, aes(x = From / 1e6, y = intact.pcnt * 100, color=Chr), lwd=0.001, position = "dodge") +
  ylab("Intact/Total LTR (%)") +
  theme(axis.title.y = element_text(size=9), legend.position = "none", axis.title.x=element_blank())
NAM.chr1_gmap_p <- ggplot() + 
  geom_density(data=NAM.rmap, aes(x=pos/1e6, , group=panel, color=panel), adjust=0.05, lwd=0.05, trim=T) + 
  xlab("B73 Chr1 position (Mb)") + ylab("Crossover Density") +
  theme(axis.title.y = element_text(size=9), legend.position = "none")

NAM.chr1_p = plot_grid(
  NAM_intact_p,
  NAM_frag_p,
  NAM.chr1_total_p,
  NAM.chr1_gypsy_p,
  NAM.chr1_copia_p,
  NAM.chr1_intact_adj_p,
  NAM.chr1_gmap_p,
  labels=c("", "", ""), ncol = 1, align = "v")

# save it to a pdf
pdf("NAM_LTR_distribution.pdf",width=6,height=10,pointsize=12, paper='special')
  NAM.chr1_p
dev.off()
