# heatmap pav

fmt_matrix <-read.csv("~/Desktop/Pan_genome_follow_up/2021_March_Pan_matrix_update/matrix_for_PAV_heatmap.csv",header = TRUE,stringsAsFactors=FALSE)

pav <- fmt_matrix[5:30]

# convert gene pav into presence, absence 
pav[pav != "NA" ] <- "Presence"
pav[is.na(pav)] <- "Absence"


# count the number of presence and sort the dataframe 
pav$count <- rowSums(pav[1:26] == "Presence")
# sort the dataset by the number of NA, decending 
library(dplyr)
sort_pav <- pav[with(pav, order(-pav$count)), ]
sequence <- seq(from =1, to =103033, by =1)
# create dataframe so tidyr can transform the dataset for heatmap 
heatmap_dataset <- cbind(sequence,sort_pav[1:26])
# reshape the dataframe and order column by genome 
heatmap_dataset_for_plot <- subset(heatmap_dataset, select = c("sequence", "Tzi8", "NC358", "NC350","Ki11","Ki3","CML333", "CML322","CML277", "CML247", "CML228","CML103","CML69","CML52","Il14H","P39","HP301","Tx303","Mo18W","M37W","Oh7B","Oh43","Ms71","M162W","Ky21","B97","B73"))

transform_hp <- melt(heatmap_dataset_for_plot, id.var = 'sequence')
#pdf("pav_heatmap.pdf", width=6,height=6,pointsize=12)

color_fill <- c("#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#FF4500","#FF4500","#DA70D6","#787878","#787878","#787878","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#FFC125")
PAV_heatmap <- ggplot(transform_hp, aes(sequence, variable)) + geom_tile(aes(fill = value)) + 
  scale_fill_manual(name="Pan-Genes (n=103033) Presence or Absence",values=c("#316879","#ff9a8d"))  + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                                                              axis.text.x=element_blank(),
                                                                                                                              axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(angle = 0, hjust = 1, colour = color_fill),
        axis.ticks.y=element_blank())  +
  theme(legend.position="bottom") 


PAV_heatmap + theme(text = element_text(size = 12)) 

