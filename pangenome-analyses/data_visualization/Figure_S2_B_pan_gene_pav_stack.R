# pan-gene type per genome and heatmap pav 

library(gdata)

plot_gene_type <- read.csv("summary_gene_type_fmt.csv")
reshaped_datafram <- pivot_wider(plot_gene_type, names_from = NAM, values_from = Gene_Count)

#plot_type_fmt <- read.csv("reshaped_gene_type_fmt.csv",header = TRUE)
plot_type_dataset_for_plot <- subset(reshaped_datafram, select = c(Gene_Type,B73,B97,Ky21,M162W,Ms71,Oh43,Oh7B,M37W,Mo18W,Tx303,HP301,P39,Il14H,CML52,CML69,CML103,CML228,CML247,CML277,CML322,CML333,Ki3,Ki11,NC350,NC358,Tzi8))
# transpose this table for supplimental table 3 
plot_type_dataset_tableS3 <- as.data.frame(t(as.matrix(plot_type_dataset_for_plot))) %>% write.csv("TableS3_pan_gene_count_per_genome.csv")

plot_type_dataset_for_plot <- trim(plot_type_dataset_for_plot, recode.factor=FALSE)


# convert matrix back to three columns for stacked plot 
library(reshape2)
plot_type_dataset_for_plot_reshape <- setNames(melt(plot_type_dataset_for_plot), c('Gene_Type', 'NAM', 'Gene_Count'))
color_fill <- c("#FFC125","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#4169E1","#787878","#787878","#787878","#DA70D6","#FF4500","#FF4500","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32","#32CD32")
levels(factor(plot_type_dataset_for_plot_reshape$Gene_Type))

gene_composition <- ggplot(plot_type_dataset_for_plot_reshape, aes(y=Gene_Count, x=factor(NAM), fill=factor(Gene_Type,levels=c("Core Gene","Near-Core Gene", "Dispensable Gene","Private Gene" )))) +
  geom_bar(position = "fill", stat = "identity") +
  ylab("Pan-Gene Composition in NAM Genome") +
  #scale_fill_npg(name="Pan Gene Type") + 
  scale_fill_manual(values = c("#DC0000FF","#3C5488FF","#4DBBD5FF","#00A087FF"),labels = c("Core Gene", "Near-Core Gene", "Dispensable Gene","Private Gene")) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab("NAM Genomes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90, hjust = 1, colour = color_fill),legend.title = element_blank()) + theme(legend.position="bottom") 


gene_composition + theme(text = element_text(size = 12))


# to get all the percentage number 
plot_type_dataset_for_plot_reshape
summ <- ddply(plot_type_dataset_for_plot_reshape, .(NAM, Gene_Type), summarize, Sum_gene=sum(Gene_Count))
gene_type_percentage <- ddply(summ, .(NAM), mutate, gene_pct = Sum_gene / sum(Sum_gene) * 100)

core_gene = gene_type_percentage %>% filter(Gene_Type == "Core Gene")
core_mean <- mean(core_gene$gene_pct)
std.error(core_gene$gene_pct)
sum(core_gene$Sum_gene)


near_core = gene_type_percentage %>% filter(Gene_Type == "Near-Core Gene")
near_core_mean <- mean(near_core$gene_pct)
std.error(near_core$gene_pct)


dispensable = gene_type_percentage %>% filter(Gene_Type == "Dispensable Gene")
mean(dispensable$gene_pct)
std.error(dispensable$gene_pct)

private = gene_type_percentage %>% filter(Gene_Type == "Private Gene")
mean(private$gene_pct)
std.error(private$gene_pct)
