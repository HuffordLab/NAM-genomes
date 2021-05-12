# pan gene matrix visulization 
library(dplyr)
library(ggsci)

pan_matrix <- read.csv("pan_gene_matrix_march_v3_all_info.csv")
subgenome_plotting <- pan_matrix[,31:33]
subgenome_plotting$Subgenome <- as.character(subgenome_plotting$Subgenome)
subgenome_plotting$Subgenome[is.na(subgenome_plotting$Subgenome)] <- "Non-Syntenic"

summary_count <-
  subgenome_plotting %>%
  count(class, Subgenome, number_genome_presence,sort = TRUE)
as.data.frame(summary_count)

gene_frequency_plot <- ggplot(as.data.frame(summary_count), aes(x = number_genome_presence, y = n)) + 
  geom_bar(aes(fill = class), position = "stack", stat = "identity",fill = "black") +
  geom_bar(aes(alpha = Subgenome), stat = "identity", fill = "grey") +
  scale_alpha_manual(values = c(0.9, 0.5, 0)) + 
  labs(x = "Number of Genomes",
       y = "Number of Pan Genes",
       fill = "Pan Gene Type",
       alpha = "Maize Subgenome") +
  scale_x_continuous(breaks = seq(1, 26, 1)) +
  ylim(0, 30000) + 
  theme(text = element_text(size = 14),legend.key = element_rect(fill = "black")) + theme_classic()

gene_frequency_plot + theme(legend.position = "none") 


####
# Piechart section 
pan_gene_pct  <-
  pan_matrix %>% select(class) %>%
  count(class,sort = TRUE)

pan_freq <- as.data.frame(pan_gene_pct)

pan_freq$class <- factor(pan_freq$class, levels = c("Core Gene", "Dispensable Gene","Private Gene","Near-Core Gene"))
pan_freq %>% mutate(percent = n/sum(n))
# class     n    percent
# 1 Dispensable Gene 51093 0.49588967
# 2        Core Gene 27910 0.27088409
# 3     Private Gene 19888 0.19302554
# 4   Near-Core Gene  4142 0.04020071

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

ggplot(pan_freq, aes(x="", y=n, fill=class))+
  geom_bar(width = 1, stat = "identity")  +scale_fill_npg() + 
  coord_polar("y") + blank_theme + theme(legend.position ="none") 

bp<- ggplotGrob(ggplot(pan_freq, aes(x="", y=n, fill=class))+
                  geom_bar(width = 1, stat = "identity")  +scale_fill_npg() + 
                  coord_polar("y") + blank_theme + theme(legend.position ="none"))



#pdf("~/Desktop/NAM/gene_frequency_pie1.pdf", width=11,height=6,pointsize=12, paper='special')
pie_chart_gene_frequency <- pie + blank_theme 


pan_gene_frequency_anchor <- gene_frequency_plot + 
  annotation_custom(
    grob = bp,
    xmin = 6,
    xmax = 20,
    ymin = 10000,
    ymax = 25000
  ) 

pan_gene_frequency_piechart <- 
  pan_gene_frequency_anchor + annotate("text", x = 8, y = 20000, label = "Core Genes: 27.09%",size = 4) + annotate("text", x = 15, y = 24000, label = "Near-Core Genes: 4.02%",size = 4) +
  annotate("text", x = 14, y = 15000, label = "Dispensable Genes: 49.59%",size =4) + annotate("text", x = 19, y = 20000, label = "Private Genes: 19.30%",size = 4) + 
  stat_summary(fun.y = sum, aes(label = ..y.., group = number_genome_presence), geom = "text",vjust=0.5, size=3.5, angle = 90,hjust =-0.1 )

pan_gene_frequency_piechart + theme(legend.position = "none") + theme(axis.title.y = element_text(size=14),
                                                                      axis.title.x = element_text(size=14),
                                                                      axis.text.y = element_text(size=14), 
                                                                      axis.text.x = element_text(size=14,angle = 90),
                                                                      legend.text = element_text(size=14),
                                                                      legend.title=element_blank())



# making the color bar 
summary_count$class <- factor(summary_count$class, levels = c("Core Gene", "Dispensable Gene","Private Gene","Near-Core Gene"))


gene_frequency_plot <- ggplot(summary_count, aes(x=number_genome_presence, y=n, fill=class)) +
  geom_bar(stat="identity")+
  
  scale_x_continuous(name ="Number of Genomes", breaks=seq(1,26,1)) +
  ylim(0,30000) + theme(text = element_text(size=20)) + 
  ylab("Number of Pan Genes") + scale_fill_npg() + theme_classic() + labs(fill = "Pan Gene Type") + theme(legend.position ="none") + 
  theme(legend.position = "none") + theme(axis.title.y = element_text(size=14),
                                          axis.title.x = element_text(size=14),
                                          axis.text.y = element_text(size=14), 
                                          axis.text.x = element_text(size=14,angle = 90),
                                          legend.text = element_text(size=14),
                                          legend.title=element_blank())



gene_frequency_plot + theme(legend.position = "none") 


