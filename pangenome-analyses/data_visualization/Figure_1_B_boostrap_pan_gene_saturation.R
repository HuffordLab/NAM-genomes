library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
# prepare dataset 
# bootstrap_stepwise
color_levels = brewer.pal(8,'Greys')

pan_gene_bs_26 = read.table('26_genome_1000bs.txt', header=F)
pan_gene_bs_26 = data.frame(count=as.vector(t(pan_gene_bs_26)),Genome=rep(1:26, 1000))
pan_gene_bs_26$Genome = as.factor(pan_gene_bs_26$Genome)
pan_gene_bs_26$group <- "26 NAM"

pan_gene_bs_tropical = read.table('tropical_1000bs.txt', header=F)
pan_gene_bs_tropical = data.frame(count=as.vector(t(pan_gene_bs_tropical)),Genome=rep(1:13, 1000))
pan_gene_bs_tropical$Genome = as.factor(pan_gene_bs_tropical$Genome)
pan_gene_bs_tropical$group <- "Tropical"

pan_gene_bs_temperate = read.table('temperate_1000bs.txt', header=F)
pan_gene_bs_temperate = data.frame(count=as.vector(t(pan_gene_bs_temperate)),Genome=rep(1:10, 1000))
pan_gene_bs_temperate$Genome = as.factor(pan_gene_bs_temperate$Genome)
pan_gene_bs_temperate$group <- "Temperate"

# all three 
joint_all <- rbind(pan_gene_bs_26,pan_gene_bs_tropical,pan_gene_bs_temperate) 
joint_all$Genome = as.factor(joint_all$Genome)

Total_pan_genes = 103033

# Fill violin with each group color 
joint_all_Fig <- ggplot(joint_all, aes(x=Genome, y=count,color=group,fill=group)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf,color = NA, ymax = Total_pan_genes*0.8),
            fill = color_levels[2]) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.8,color = NA, ymax = Total_pan_genes*0.9),
            fill = color_levels[3]) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.9,color = NA, ymax = Total_pan_genes*0.95),
            fill = color_levels[4]) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.95,color = NA, ymax = Inf),
            fill = color_levels[6]) +
  geom_violin() + #this layer on the top
  scale_colour_manual(values = c( "#000000","#4169E1","#32CD32"),labels = c("26 NAM","Temperate", "Tropical")) +
  scale_fill_manual(values = c( "#000000","#4169E1","#32CD32"),labels = c("26 NAM","Temperate", "Tropical")) +
  scale_y_continuous(limits=c(40000, 105000)) +
  labs(x ="Number of Genomes", y = "Number of Pan-Genes")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14,angle = 90),
        legend.text = element_text(size=14),
        legend.title=element_blank())

require(scales)
joint_all_Fig + scale_y_continuous(labels = comma, sec.axis = sec_axis(~ . *100/Total_pan_genes, breaks = c(50,60,70,80,85,90,95,100),name="Percentage of Pan-Genes (%)")) + 
  theme(strip.text.x = element_blank(), 
        strip.background = element_rect(colour=color_levels[2], fill=color_levels[2]),
        legend.position=c(.8,.3)
  )
