# PAV heatmap 

setwd("~/Desktop/NAM_PAN_GENOME/pan_genome_nov 2/QC_set/pan_gene_frequency_within_genome_pav/")
pav <- read.csv('pav_heatmap_matrix_fmt_final.csv')
dim(pav)
#after reformatting
# count the number of presence and sort the dataframe 
pav$count <- rowSums(pav[1:26] == "Presence")
# sort the dataset by the number of NA, decending 
library(dplyr)
sort_pav <- pav[with(pav, order(-pav$count)), ]
write.csv(sort_pav[1:27], file = "pav_count.csv")
dim(sort_pav)
sequence <- seq(from =1, to =103538, by =1)
# create dataframe so tidyr can transform the dataset for heatmap 
heatmap_dataset <- cbind(sequence,sort_pav[1:26])
# reshape the dataframe and order column by genome 
heatmap_dataset_for_plot <- subset(heatmap_dataset, select = c("sequence", "Tzi8", "NC358", "NC350","Ki11","Ki3","CML333", "CML322","CML277", "CML247", "CML228","CML103","CML69","CML52","Il14H","P39","HP301","Tx303","Mo18W","M37W","Oh7B","Oh43","Ms71","M162W","Ky21","B97","B73"))

transform_hp <- melt(heatmap_dataset_for_plot, id.var = 'sequence')
#pdf("pav_heatmap.pdf", width=6,height=6,pointsize=12)

color_fill <- c("#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#FF4500","#FF4500","#DA70D6","#787878","#787878","#787878","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#FFC125")
PAV_heatmap <- ggplot(transform_hp, aes(sequence, variable)) + geom_tile(aes(fill = value)) + 
  scale_fill_manual(name="Pan-Genes (n=103538) Presence or Absence",values=c("#316879","#ff9a8d"))  + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                                                              axis.text.x=element_blank(),
                                                                                                                              axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(angle = 0, hjust = 1, colour = color_fill),
        axis.ticks.y=element_blank())  +
  theme(legend.position="bottom") 


PAV_heatmap + theme(text = element_text(size = 12)) 
dev.off()
