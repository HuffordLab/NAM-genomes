# Phylostrata analysis 

library(plotrix)
library(tidyr)
library(matrixStats)
library(data.table)
library(dplyr)
library(ggplot2)
library(lemon)


#tandem copy number matrix 
tandem_copy_list <- read.csv("pan_gene_copy_number.txt",sep = "\t",header=FALSE)
colnames(tandem_copy_list) <- c("Query_gene","NAM_genome","Tandem")
#reshape the dataset 
NAM_base_copy_number <- subset(tandem_copy_list,tandem_copy_list$NAM_genome != "pan_gene" & tandem_copy_list$NAM_genome !="Query_gene" )
#fill_info[rowSums(tandem_copy_list=="")!=ncol(tandem_copy_list), ]
reshaped_datafram_copy_number <- pivot_wider(NAM_base_copy_number, names_from = NAM_genome, values_from = Tandem)

# reorder column name before joint the dataset 
matrix_for_copy_count= subset(reshaped_datafram_copy_number, select = c(Query_gene,B73,Tzi8,Ky21,M162W,Ms71,Oh7B,Oh43,M37W,Mo18W,NC350,HP301,Il14H,P39,CML52,CML69,Ki11,CML228,CML247,CML277,CML322,CML333,Ki3,CML103,Tx303,NC358,B97))
# add in pan gene type information, 0 will be considered as NA, missing gene
matrix_for_copy_count$occurance_genome <- rowSums(matrix_for_copy_count[2:27]!=0)

#classify geneome catogary 
matrix_for_copy_count$class<-ifelse(matrix_for_copy_count$occurance_genome ==26,rr2<-"Core Gene",
                                    ifelse(matrix_for_copy_count$occurance_genome >23 &matrix_for_copy_count$occurance_genome<26,rr2<-"Near-Core Gene",
                                           ifelse(matrix_for_copy_count$occurance_genome>1 &matrix_for_copy_count$occurance_genome<24,rr2<-"Dispensable Gene",
                                                  rr2<-"Private Gene")))

# write as numeric pan matrix 
write.csv(matrix_for_copy_count, "numeric_pan_matrix.csv")

# sum of the copy number per row
copy_sum <- matrix_for_copy_count %>%
  mutate(Total = select(., 2:27) %>% rowSums(na.rm = TRUE)) 

summary(copy_sum$Total)

# adding in phylostrata info 
base_pan_phylo <- read.csv("phylo_info.txt",sep = "\t",header =FALSE)
names(base_pan_phylo) <- c("Query_gene","Evi_ab","Phylo")

# create join matrix with copy number, pan gene type, and phylostrata information

join_matrix <- left_join(copy_sum,base_pan_phylo) #%>% left_join(pan_gene_canonical)

# subset the last four columns for the phylostrata figure visulization 

keys <- colnames(join_matrix[29:32])[!grepl('Total',colnames(join_matrix[29:32]))]
X <- as.data.table(join_matrix[29:32])
sum_total <- X[,list(Total_present= sum(Total)),keys]

sum(sum_total$Total_present)
# [1] 1383856 genes included here

summary(sum_total$Evi_ab)
aggregate(sum_total$Total_present, by=list(Evi_ab=sum_total$Evi_ab), FUN=sum)


sum_total$class <- factor(sum_total$class, levels=c("Core Gene","Near-Core Gene","Dispensable Gene","Private Gene"))
head(sum_total)

color_fill <- c("darkred", "darkblue", "gold", "darkgreen")
a1 <- col2rgb(color_fill)
a2 <- rgb2hsv(a1)

stack_all <- sum_total %>% 
  ggplot(aes(y=Total_present, x=factor(class), fill=factor(Phylo,levels=c("Viridiplanteae","Poaceae", "Andropogoneae","Maize" )))) +
  geom_bar(position = "fill", stat = "identity") + theme(legend.title = element_blank()) + ggtitle("Full") + 
  scale_fill_manual(values = hsv(a2[1,], a2[2,]*0.4, a2[3,]),labels = c("Viridiplanteae", "Poaceae", "Andropogoneae","Maize"))+ 
  theme_classic()+ ylab("Propotion in Phylostrata") + 
  theme(axis.title.x=element_blank(),axis.text.x = element_text(color="black",size=12, angle=90)) + labs(fill="") 

# abinitio 
stack_ab <- sum_total %>% filter(Evi_ab =="abinitio") %>%
  ggplot(aes(y=Total_present, x=factor(class), fill=factor(Phylo,levels=c("Viridiplanteae","Poaceae", "Andropogoneae","Maize" )))) +
  scale_fill_manual(values = hsv(a2[1,], a2[2,]*0.4, a2[3,]),labels = c("Viridiplanteae", "Poaceae", "Andropogoneae","Maize"))+ 
  geom_bar(position = "fill", stat = "identity") + theme(legend.title = element_blank()) + ggtitle("Ab initio") + theme_classic()+ theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank(),axis.text.x = element_text(color="black",size=12, angle=90)) + labs(fill="")


stack_evi <- sum_total %>% filter(Evi_ab =="evidence") %>%
  ggplot(aes(y=Total_present, x=factor(class), fill=factor(Phylo,levels=c("Viridiplanteae","Poaceae", "Andropogoneae","Maize" )))) +
  scale_fill_manual(values = hsv(a2[1,], a2[2,]*0.4, a2[3,]),labels = c("Viridiplanteae", "Poaceae", "Andropogoneae","Maize"))+ 
  geom_bar(position = "fill", stat = "identity") + theme(legend.title = element_blank()) + ggtitle("Evidence") + theme_classic()+ theme(axis.title.y=element_blank()) + 
  theme(axis.title.x=element_blank(),axis.text.x = element_text(color="black",size=12, angle=90)) + labs(fill="")

grid_arrange_shared_legend(stack_all, stack_ab,stack_evi, nrow = 1) 
