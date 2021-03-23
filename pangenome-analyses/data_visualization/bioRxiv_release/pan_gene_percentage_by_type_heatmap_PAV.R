library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(ggsci)

setwd("~/Desktop/pan_genome_nov 2/QC_set/pan_gene_frequency_within_genome_pav/")

pan_matrix <- read.csv("final_pan_matrix_for_visualization.csv")
# count how many NA per pan genes have
pan_matrix$na_count <- apply(pan_matrix, 1, function(x) sum(is.na(x)))

# plot bar pan gene by NAM line 
# making types of pan genes  
pan_matrix$pan_gene <-ifelse(pan_matrix$na_count<1,rr2<-"Core Gene",
                             ifelse(pan_matrix$na_count>0&pan_matrix$na_count<3,rr2<-"Soft Gene",
                                    ifelse(pan_matrix$na_count>2&pan_matrix$na_count<25,rr2<-"Dispensible Gene",
                                           rr2<-"Private Gene")))


# store dataframe that has the gene class information 
write.csv(pan_matrix, file = "pav_matrix_with_NA_count.csv")

#reformat the csv in excel, replace query gene by the pan_gene type 
bar_plot <- read.csv("pav_matrix_with_NA_count_final.csv",header = TRUE,stringsAsFactors=FALSE)
#remove gene name and replace by pan_gene_type
pan_gene_type_matrix<- sapply(bar_plot[-1], function(x) {ind <- which(x!="NA"); x[ind] = bar_plot[ind,2]; return(x)})
write.csv(pan_gene_type_matrix,file="pav_heatmap_matrix.csv")
summary_gene_type <- as.data.frame(summary(pan_gene_type_matrix))
write.csv(summary_gene_type, file="summary_gene_type.csv")

#format the summary file, give header and remove lines with NA count (a total of 104 lines left)
# this generate the stacked bar plot with % of genes in the genome instead of solo gene number count 
library(tidyr)

plot_gene_type <- read.csv("summary_gene_type_fmt_final.csv")
reshaped_datafram <- pivot_wider(plot_gene_type, names_from = NAM, values_from = Gene_Count)

#plot_type_fmt <- read.csv("reshaped_gene_type_fmt.csv",header = TRUE)
plot_type_dataset_for_plot <- subset(reshaped_datafram, select = c(Gene_Type,B73,B97,Ky21,M162W,Ms71,Oh43,Oh7B,M37W,Mo18W,Tx303,HP301,P39,Il14H,CML52,CML69,CML103,CML228,CML247,CML277,CML322,CML333,Ki3,Ki11,NC350,NC358,Tzi8))

# convert matrix back to three columns for stacked plot 
library(reshape2)
plot_type_dataset_for_plot_reshape <- setNames(melt(plot_type_dataset_for_plot), c('Gene_Type', 'NAM', 'Gene_Count'))
color_fill <- c("#FFC125","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#787878","#787878","#787878","#DA70D6","#FF4500","#FF4500","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32")
levels(factor(plot_type_dataset_for_plot_reshape$Gene_Type))

gene_composition <- ggplot(plot_type_dataset_for_plot_reshape, aes(y=Gene_Count, x=factor(NAM), fill=factor(Gene_Type,levels=c("Core Gene       ","Soft Gene       ", "Dispensible Gene","Private Gene    " )))) +
  geom_bar(position = "fill", stat = "identity") +
  #geom_text(aes(label = Gene_Count), position = position_fill(vjust = 0.5)) +
  ylab("Percentage of Genes in the Pan Genome") +
  #scale_fill_npg(name="Pan Gene Type") + 
  scale_fill_manual(values = c("#DC0000FF","#3C5488FF","#4DBBD5FF","#00A087FF"),labels = c("Core Gene", "Near-Core Gene", "Dispensable Gene","Private Gene")) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab("NAM Genomes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45, hjust = 1, colour = color_fill),legend.title = element_blank()) + theme(legend.position="bottom") 


gene_composition + theme(text = element_text(size = 12))



# repalce pan gene type with presence/absence before moving into the heatmap section 

sed 's/Dispensible Gene/Presence/g; s/bc/ab/g; s/Soft Gene/Presence/g; s/Core Gene/Presence/g; s/~~/bc/g; s/NA/Absence/g; s/Private Gene/Presence/g' pav_heatmap_matrix.csv > pav_heatmap_matrix_fmt_final.csv


# PAV heatmap 
pav <- read.csv('pan_gene_frequency_within_genome_pav/pav_heatmap_matrix_fmt_final.csv')
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
  scale_fill_manual(name="Pan Genes (n=103538) Presence or Absence",values=c("#316879","#ff9a8d"))  + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                                                              axis.text.x=element_blank(),
                                                                                                                              axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(angle = 45, hjust = 1, colour = color_fill),
        axis.ticks.y=element_blank())  +
  theme(legend.position="bottom") 


PAV_heatmap + theme(text = element_text(size = 12)) 
dev.off()

