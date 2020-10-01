library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(ggsci)

setwd("~/Desktop/NAM_pan_genome_final/909090/pan_gene_frequency/")
# read in collasped final matrix, add the pan gene as the firs column "pan_gene_##" 
pan_gene_matrix <- read.csv("final_pan_matrix_rename.csv",header = TRUE,stringsAsFactors=FALSE)
# rename the matrix, for NAM lines with the pan gene present, the NAM gene will be replaced by the pan gene ID. If the gene is missing, it remains to be NA
pan_gene_matrix_rename_by_pan_ID <- sapply(pan_gene_matrix[,-1], function(x) {ind <- which(x!="NA"); x[ind] = pan_gene_matrix[ind,1]; return(x)})

# count how many NA per pan genes have
pan_gene_matrix$na_count <- apply(pan_gene_matrix, 1, function(x) sum(is.na(x)))
# get the frequency table based on the NA count 
pan_gene_frequency <- as.matrix(table(pan_gene_matrix$na_count))
write.csv(pan_gene_frequency, file = "pan_gene_frequency.csv") # this dataset includes NA frequency, which helps to visualize the pan gene frequency by class we define


gene_freq <- read.csv("pan_gene_frequency.csv", header = TRUE)
colnames(gene_freq) <- c("Number_of_Genomes", "Pan_gene")
#classify geneome catogary 
gene_freq$class<-ifelse(gene_freq$Number_of_Genomes ==0,rr2<-"Core Gene",
                        ifelse(gene_freq$Number_of_Genomes>0 &gene_freq$Number_of_Genomes<3,rr2<-"Softcore Gene",
                               ifelse(gene_freq$Number_of_Genomes>2 &gene_freq$Number_of_Genomes<25,rr2<-"Dispensible Gene",
                                      rr2<-"Private Gene")))
# add a column with reversed genomue value for visualization 
genome_presence <- rev(gene_freq$Number_of_Genomes+1)


#rename the figure 
gene_freq_for_visualization <- cbind(gene_freq,genome_presence)
write.csv(gene_freq_for_visualization,file="gene_freq_for_visualization.csv")
# change column name into Number.of.Genomes and Pan_gene
# pdf("~/Desktop/NAM/pan_gene_frequency1.pdf", width=11,height=6,pointsize=12, paper='special')

gene_frequency_plot <- ggplot(gene_freq_for_visualization, aes(x=genome_presence, y=Pan_gene, fill=class)) +
  geom_bar(stat="identity")+
  
  scale_x_continuous(name ="Number of Genomes", breaks=seq(0,26,1)) +
  ylim(0,30000) + theme(text = element_text(size=20)) + 
  ylab("Number of Pan Genes") + scale_fill_npg() + theme_classic() + labs(fill = "Pan Gene Type") + theme(legend.position ="none")


gene_frequency_plot + theme(text = element_text(size = 12))
### add the piechart in the exisiting plot 
# piechart version 
blank_theme <- theme_minimal() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title=element_text(size=0, face="bold")
  ) 

bp<- ggplotGrob(ggplot(gene_freq_for_visualization, aes(x="genome_presence", y=Pan_gene, fill=class))+
                  geom_bar(width = 1, stat = "identity")  +scale_fill_npg() + 
                  coord_polar("y") + blank_theme + labs(fill = "Pan Gene Type")+ theme(legend.position ="none") )


#pdf("~/Desktop/NAM/gene_frequency_pie1.pdf", width=11,height=6,pointsize=12, paper='special')
pie_chart_gene_frequency <- pie + blank_theme 


pan_gene_frequency_anchor <- gene_frequency_plot + 
  annotation_custom(
    grob = bp,
    xmin = 2,
    xmax = 25,
    ymin = 5000,
    ymax = 27000
  ) 

pan_gene_frequency_piechart <- pan_gene_frequency_anchor + annotate("text", x = 8, y = 18000, label = "Core Genes: 27.13%",size = 4.5) + annotate("text", x = 16, y = 26000, label = "Softcore Genes: 4.00%",size = 4.5) +
  annotate("text", x = 14, y = 12000, label = "Dispensible Genes: 49.08%",size =4.5) + annotate("text", x = 20, y = 18000, label = "Priviate Genes: 19.79%",size = 4.5) +  theme(text = element_text(size = 12)) + geom_text(aes(label=Pan_gene), vjust=0.5, size=4, angle = 90,hjust =-0.1) 


pdf("pan_gene_frequency_pdf_print.pdf", width=7,height=6,pointsize=12)
pan_gene_frequency_piechart
dev.off()
