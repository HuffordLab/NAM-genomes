library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(ggsci)

setwd("~/Desktop/NAM_PAN_GENOME/pan_genome_nov 2/QC_set/pan_gene_frequency_within_genome_pav/")

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
  ylab("Pan-Gene Composition in NAM Genome") +
  #scale_fill_npg(name="Pan Gene Type") + 
  scale_fill_manual(values = c("#DC0000FF","#3C5488FF","#4DBBD5FF","#00A087FF"),labels = c("Core Gene", "Near-Core Gene", "Dispensable Gene","Private Gene")) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab("NAM Genomes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90, hjust = 1, colour = color_fill),legend.title = element_blank()) + theme(legend.position="bottom") 


gene_composition + theme(text = element_text(size = 12))
