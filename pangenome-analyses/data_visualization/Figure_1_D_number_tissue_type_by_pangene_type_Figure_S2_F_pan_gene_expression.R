library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(lemon)
library(ggpubr)
library(ggsignif)

Tissue_8_expression_matrix <- read.csv(file = "expression_26_NAM_8_tissues_rm_header.csv", header = FALSE)

# adding column name 
colnames(Tissue_8_expression_matrix) <- c("PanGeneID","Geneid","8DAS_root","8DAS_shoot","R1_anther","V11_base",	"V11_middle","V11_tip",	"V18_ear","V18_tassel","number_of_tissue_expression", "pan-gene class")		

# adding mean expression value per pan gene. When pan gene class is NA, meaning in the updated pan gene matrix no longer have that gene (merged via the recent pan matrix update)
Tissue_8_expression_matrix_clean = Tissue_8_expression_matrix %>% filter(`pan-gene class` != "NA")

Tissue_8_expression_matrix_clean$mean_expression <- rowMeans(Tissue_8_expression_matrix_clean[3:10],na.rm = TRUE)
table(Tissue_8_expression_matrix_clean$`pan-gene class`)
# Core Gene Dispensable Gene   Near-Core Gene     Private Gene 
# 725489           393801           102195            19888 


##########  Pan gene expression RPKM    #############
# for expression RPKM, all pan gene that has expression level NaN is excluded
RPKM_visluization = Tissue_8_expression_matrix_clean %>% filter(mean_expression != "NaN")

table(RPKM_visluization$`pan-gene class`)
#       Core Gene   Near-Core Gene Dispensable Gene     Private Gene 
#       565016            46015            70741             6376 
# Non-Log transform 
RPKM_visluization$`pan-gene class` <- factor(RPKM_visluization$`pan-gene class`,levels = c("Core Gene","Near-Core Gene","Dispensable Gene","Private Gene"))


pan_expression_RPKM<- RPKM_visluization %>% filter(mean_expression <25) %>% 
  ggplot(aes(x=`pan-gene class`,y=mean_expression,fill=`pan-gene class`),palette = "jco") + geom_boxplot(outlier.shape = NA) + theme_classic() + 
  scale_fill_manual(values = c("#DC0000FF","#3C5488FF","#4DBBD5FF","#00A087FF"),labels = c("Core Gene", "Near-Core Gene", "Dispensable Gene","Private Gene")) + 
  labs(x = "Pan-Gene Type", y = "RPKM",fill = "Pan-Gene Type") +   theme(strip.text.x = element_blank(), 
                                                                         strip.background = element_rect(colour=color_levels[2], fill=color_levels[2]),
                                                                         legend.position=c(.8,0.90)) + 
  geom_signif(
    comparisons = list(c("Core Gene", "Near-Core Gene")),
    map_signif_level = FALSE, textsize = 4, y_position = 21,
    test="t.test", test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE)) + 
  geom_signif(
    comparisons = list(c("Near-Core Gene", "Dispensable Gene")),
    map_signif_level = FALSE, textsize = 4,test="t.test", test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE),
    y_position = 17
  ) + geom_signif(
    comparisons = list(c("Dispensable Gene", "Private Gene")),
    map_signif_level = FALSE, textsize = 4,test="t.test", test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE),
    y_position = 13
  ) + geom_signif(
    comparisons = list(c("Near-Core Gene", "Private Gene")),
    map_signif_level = FALSE, textsize = 4,test="t.test", test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE),
    y_position = 14.5) + 
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=0),
        legend.text = element_text(size=14),
        legend.title=element_blank()) 




#### number of tissue type per pan gene ##############

# for this analysis, all the pan gene will be included. no expression will be counted as absent in the tissue type 

tissue_type_expression  <- Tissue_8_expression_matrix_clean %>% select(`pan-gene class`,number_of_tissue_expression)

table(tissue_type_expression$`pan-gene class`)

tissue_type_expression_melt = melt(table(tissue_type_expression)) 
tissue_type_expression_melt$number_of_tissue_expression <- factor(tissue_type_expression_melt$number_of_tissue_expression)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=12)
  )

# core 
Core_pie <- tissue_type_expression_melt %>% filter(`pan-gene class` == "Core Gene") %>% 
  ggplot(aes(x="", y=value, fill=number_of_tissue_expression)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + scale_fill_grey(start=0.8, end=0.2,name = "Number of Tissue Types") + blank_theme + theme(axis.text.x=element_blank()) + 
  ggtitle("Core Gene") + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.title=element_text(size=14),
        legend.text = element_text(size=14),plot.title = element_text(hjust = 0.5)) 

# (n=725,489)


Near_Core_pie <- tissue_type_expression_melt %>% filter(`pan-gene class` == "Near-Core Gene") %>% 
  ggplot(aes(x="", y=value, fill=number_of_tissue_expression)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + scale_fill_grey(start=0.8, end=0.2) + blank_theme + theme(axis.text.x=element_blank()) + 
  ggtitle("Near-Core Gene") + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.title=element_text(size=14),
        legend.text = element_text(size=14),plot.title = element_text(hjust = 0.5)) 

# (n=102,195)

Dispensable_pie <- tissue_type_expression_melt %>% filter(`pan-gene class` == "Dispensable Gene") %>% 
  ggplot(aes(x="", y=value, fill=number_of_tissue_expression)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + scale_fill_grey(start=0.8, end=0.2) + blank_theme + theme(axis.text.x=element_blank()) +
  ggtitle("Dispensable Gene") + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.title=element_text(size=14),
        legend.text = element_text(size=14),plot.title = element_text(hjust = 0.5)) 

#  (n=393,801)


Private_pie <- tissue_type_expression_melt %>% filter(`pan-gene class` == "Private Gene") %>% 
  ggplot(aes(x="", y=value, fill=number_of_tissue_expression)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + scale_fill_grey(start=0.8, end=0.2) + blank_theme + theme(axis.text.x=element_blank()) + 
  ggtitle("Private Gene") + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.title=element_text(size=14),
        legend.text = element_text(size=14),plot.title = element_text(hjust = 0.5)) 

#  (n=19,888)

grid_arrange_shared_legend(Core_pie,Near_Core_pie,Dispensable_pie,Private_pie,nrow=1) 


#######blow numbers are used for response to reviewers
# core gene in 8 tissue types 
# get percentage of core gene in 8 tissue type 
total_core = tissue_type_expression %>% filter(`pan-gene class`=="Core Gene")
dim(total_core)
# 725489 
# total in core in 8 tissue
total_core_in_8 = tissue_type_expression %>% filter(`pan-gene class`=="Core Gene") %>% filter(number_of_tissue_expression== 8) 
dim(total_core_in_8)
# 287125  

# 287125/725489 = 0.3957675 

# Dispensable and private that are in at least than 1 tissue type 
total_dispensable = tissue_type_expression %>% filter(`pan-gene class`=="Dispensable Gene")
dim(total_dispensable)
# 393801 
# total in core in 8 tissue
total_dispensable_in_8 = tissue_type_expression %>% filter(`pan-gene class`=="Dispensable Gene") %>% filter(number_of_tissue_expression >0 ) 
dim(total_dispensable_in_8)
# 70741  
# 70741/393801 = 0.1796364




# Private and private that are in at least than 1 tissue type 
total_private = tissue_type_expression %>% filter(`pan-gene class`=="Private Gene")
dim(total_private)
# 19888 
# total in core in 8 tissue
total_private_in_at_least1 = tissue_type_expression %>% filter(`pan-gene class`=="Private Gene") %>% filter(number_of_tissue_expression >0 ) 
dim(total_private_in_at_least1)
# 6376  
# 6376/19888 = 0.3205953

