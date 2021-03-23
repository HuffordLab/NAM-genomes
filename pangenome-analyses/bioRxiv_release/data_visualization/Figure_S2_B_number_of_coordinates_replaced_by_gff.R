setwd("~/Desktop/NAM_PAN_GENOME/pan_genome_nov 2/QC_set/gff_intersecting/")
name_replace <- read.csv("gff_intersect_name_replace.csv")
NA_fill_for_plot <- subset(name_replace, select = c(NAM,B73,B97,Ky21,M162W,Ms71,Oh43,Oh7B,M37W,Mo18W,Tx303,HP301,P39,Il14H,CML52,CML69,CML103,CML228,CML247,CML277,CML322,CML333,Ki3,Ki11,NC350,NC358,Tzi8))
library(reshape2)
NA_fill_for_plot_reshape <- setNames(melt(NA_fill_for_plot), c("class" ,'NAM', 'gff_fill_Count'))


color_fill <- c("#FFC125","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#787878","#787878","#787878","#DA70D6","#FF4500","#FF4500","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32")


ggplot(NA_fill_for_plot_reshape, aes(x=NAM, y=gff_fill_Count)) + 
  geom_bar(stat = "identity", width=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90, hjust = 1, colour = color_fill)) + 
  xlab("NAM Genomes") + ylab("Number of Gmap Coordinates Replaced by Gene Model") + theme(text = element_text(size = 12)) 
