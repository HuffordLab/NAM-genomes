# R script to extract files with gene pan_gene ID
# change gene ID into pan gene ID (this step was done in excel)

setwd("~/Desktop/pan_genome_nov 2/QC_set/bootstrap/")

pan_figure <- read.csv("final_pan_matrix_for_visualization.csv",header = TRUE,stringsAsFactors=FALSE)
pan_figure2<- sapply(pan_figure[,-1], function(x) {ind <- which(x!="NA"); x[ind] = pan_figure[ind,1]; return(x)})



# output gene list for each name genome in R 

setwd("~/Desktop/pan_genome_nov 2/QC_set/bootstrap/ID")

col_names <- colnames(pan_figure2)
for (i in 1 : ncol(pan_figure2)) {
  nam_pan_gene_id <- as.matrix(pan_figure2[,i])
  write.csv(nam_pan_gene_id, file = paste(col_names[i], "_pan_id.csv", sep = ""))
}
write.csv(pan_figure2, file = "pan_gene_names.csv")

# commands below were done in Unix to prepare list of genes used for bootstrapping
for i in *.csv ; do
cut -d ',' -f 2 "$i"  | sed 's/"//g' | grep -v pan_all | grep -v V1 > $(basename "$i")_rename_list.txt
done 

# bootstrap method 
for k in {1..10}; do for j in {1..100}; do for i in `ls *list.txt|shuf`; do cat $i >> temp.$j.$k; sort -u temp.$j.$k|wc -l; done |bash /home/hirschc1/qiuxx221/nam_pan_genome/bootstrap/transpose.sh - > result.$j.$k; rm temp.$j.$k; done & done


# transpose.sh
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}'





# bootstrap figure visualization 
# bootstrap figure visualization 

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
setwd("~/Desktop/NAM_PAN_GENOME/pan_genome_nov 2/QC_set/bootstrap/")

# Panel A Pan Gene 
##  use pan_gene_bs_plot_pcnt for figure grabbing 

## readin bootstrap pan gene data
color_levels = brewer.pal(8,'Greys')
pan_gene_bs = read.table('Dec_bootstrap_1000.txt', header=F)
pan_gene_bs = data.frame(count=as.vector(t(pan_gene_bs)),Genome=rep(1:26, 1000))
pan_gene_bs$Genome = as.factor(pan_gene_bs$Genome)
Total_pan_genes = 103538

pan_gene_bs_plot_pcnt = ggplot(pan_gene_bs, aes(x=Genome, y=count)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Total_pan_genes*0.8),
            fill = color_levels[2]) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.8, ymax = Total_pan_genes*0.9),
            fill = color_levels[3]) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.9, ymax = Total_pan_genes*0.95),
            fill = color_levels[4]) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Total_pan_genes*0.95, ymax = Inf),
            fill = color_levels[6]) + #these layers at the bottom
  geom_violin() + #this layer on the top
  scale_y_continuous(limits=c(40000, 105000)) +
  scale_y_continuous(sec.axis = sec_axis(~ . *100/Total_pan_genes, breaks = c(50,60,70,80,85,90,95,100),name="Percentage of Pan Genes (%)")) +
  labs(x ="Number of Genomes", y = "Number of Pan Genes")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14,angle = 90),
        legend.text = element_text(size=14),
        legend.title=element_blank())

require(scales)
#pdf("pav_final.pdf", width=12,height=4,pointsize=12)
pan_gene_bs_plot_pcnt + scale_y_continuous(labels = comma, sec.axis = sec_axis(~ . *100/Total_pan_genes, breaks = c(50,60,70,80,85,90,95,100),name="Percentage of Pan Genes (%)"))
#dev.off()
