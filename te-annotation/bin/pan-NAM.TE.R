#update.packages()
library(ggplot2)
library(scales)
library(RColorBrewer)


datadir="/Users/oushujun/Google Drive/study/Maize research/NAM/TE"
setwd(datadir)


########################
#### basic settings ####
########################
## readin bootstrap pan-TE data
pan_TE_bs = read.table('pan_TE_bootstrap1000.summary26.txt', header=F)
pan_TE_bs = data.frame(count=as.vector(t(pan_TE_bs)),Genome=rep(1:26, 1000))

#percentage of TE fams to the pan-TE lib (21358 fl-TE)
mean(pan_TE_bs$count[pan_TE_bs$Genome==2])/21358
mean(pan_TE_bs$count[pan_TE_bs$Genome==5])/21358
mean(pan_TE_bs$count[pan_TE_bs$Genome==10])/21358

## plot size of TE families
pan_TE_bs$Genome = as.factor(pan_TE_bs$Genome)
pan_TE_bs_plot = ggplot(pan_TE_bs, aes(x=Genome, y=count)) + 
  geom_violin() +
  labs(x ="Number of Genomes", y = "Number of TE families") +
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.title.x = element_text(size=18, face="bold"),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_blank())
pan_TE_bs_plot

# add percentage to the second axis and change background colors
show_col(brewer.pal(8,'YlGn'))
color_levels = brewer.pal(8,'YlGn')
Total_flTE_fam = 21358
pan_TE_bs_plot_pcnt = ggplot(pan_TE_bs, aes(x=Genome, y=count)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Total_flTE_fam*0.8),
            fill = color_levels[2], alpha = 0.01) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_flTE_fam*0.8, ymax = Total_flTE_fam*0.9),
            fill = color_levels[3], alpha = 0.01) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_flTE_fam*0.9, ymax = Total_flTE_fam*0.95),
            fill = color_levels[4], alpha = 0.01) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_flTE_fam*0.95, ymax = Inf),
            fill = color_levels[5], alpha = 0.01) + #these layers at the bottom
  geom_violin() + #this layer on the top
  scale_y_continuous(sec.axis = sec_axis(~ . *100/Total_flTE_fam, breaks = c(75,80,85,90,95,100),name="Percentage of TE families (%)")) +
  labs(x ="Number of Genomes", y = "Number of TE families") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.title.x = element_text(size=18, face="bold"),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_blank())
pan_TE_bs_plot_pcnt

#output
pdf("pan-NAM.TE.pdf", width=11,height=6,pointsize=12, paper='special')
  pan_TE_bs_plot
  pan_TE_bs_plot_pcnt
dev.off()
