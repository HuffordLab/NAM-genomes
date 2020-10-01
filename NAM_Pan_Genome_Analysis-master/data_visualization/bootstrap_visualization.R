library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(ggsci)
# prepare dataset 
# bootstrap_stepwise
setwd("~/Desktop/NAM_pan_genome_final/909090/bootstrap/")

# Panel A Pan Gene 
##  use pan_gene_bs_plot_pcnt for figure grabbing 

## readin bootstrap pan gene data
color_levels = brewer.pal(8,'YlGn')
pan_gene_bs = read.table('', header=F)
pan_gene_bs = data.frame(count=as.vector(t(pan_gene_bs)),Genome=rep(1:26, 1000))
pan_gene_bs$Genome = as.factor(pan_gene_bs$Genome)
Total_pan_genes = 102660

pan_gene_bs_plot_pcnt = ggplot(pan_gene_bs, aes(x=Genome, y=count)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Total_pan_genes*0.8),
            fill = color_levels[2], alpha = 0.01) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.8, ymax = Total_pan_genes*0.9),
            fill = color_levels[3], alpha = 0.01) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.9, ymax = Total_pan_genes*0.95),
            fill = color_levels[4], alpha = 0.01) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.95, ymax = Inf),
            fill = color_levels[5], alpha = 0.01) + #these layers at the bottom
  geom_violin() + #this layer on the top
  scale_y_continuous(sec.axis = sec_axis(~ . *100/Total_pan_genes, breaks = c(50,60,70,80,85,90,95,100),name="Percentage of Pan Genes (%)")) +
  labs(x ="Number of Genomes", y = "Number of Pan Genes")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title=element_blank())

require(scales)
#pdf("pav_final.pdf", width=12,height=4,pointsize=12)
pan_gene_bs_plot_pcnt + scale_y_continuous(labels = comma, sec.axis = sec_axis(~ . *100/Total_pan_genes, breaks = c(50,60,70,80,85,90,95,100),name="Percentage of Pan Genes (%)"))
#dev.off()

