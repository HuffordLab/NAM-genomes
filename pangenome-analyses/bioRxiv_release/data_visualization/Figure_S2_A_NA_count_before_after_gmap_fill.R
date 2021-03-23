library(ggplot2)
library(reshape2)
setwd("~/Desktop/NAM_PAN_GENOME/pan_genome_nov 2/QC_set/NA_before_after/")

NA_count <- read.csv("NA_before_after_fill.csv")
plot_type_dataset_for_plot_reshape <- setNames(melt(NA_count), c('NAM', 'NA_status', 'NA_count'))
reshaped_datafram <- pivot_wider(plot_type_dataset_for_plot_reshape, names_from = NAM, values_from = NA_count)
plot_type_dataset_for_plot <- subset(reshaped_datafram, select = c(NA_status,B73,B97,Ky21,M162W,Ms71,Oh43,Oh7B,M37W,Mo18W,Tx303,HP301,P39,Il14H,CML52,CML69,CML103,CML228,CML247,CML277,CML322,CML333,Ki3,Ki11,NC350,NC358,Tzi8))
plot_type_dataset_for_plot_reshape <- setNames(melt(plot_type_dataset_for_plot), c("NA_status","NAM", "NA_count"))

plot_type_dataset_for_plot_reshape$NA_status = factor(plot_type_dataset_for_plot_reshape$NA_status, levels = c("Before_fill","Post_fill"), ordered = TRUE)


color_fill <- c("#FFC125","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#787878","#787878","#787878","#DA70D6","#FF4500","#FF4500","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32")
ggplot(plot_type_dataset_for_plot_reshape, aes(NAM, NA_count, fill=NA_status,.desc = FALSE))+ geom_bar(stat = "identity", position = 'dodge') + 
  xlab("NAM Genomes") + ylab("Number of Absent Pan-Genes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90, hjust = 1, colour = color_fill)) +
  scale_fill_manual(name="", values = c("grey","black"),labels=c("Before coordinate filling","After coordinate filling")) + 
  ylim(0,80000) +  theme(text = element_text(size = 12)) + theme(legend.position="top") 
