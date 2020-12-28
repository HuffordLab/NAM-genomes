library(reshape2)
library(viridis)
#install.packages("ggpubr")
library(ggpubr)

setwd('/Users/jianingliu/Downloads/')
nam_centromere <- data.frame(read.table('NAM.centro.coords.assembly.status.content.sum'))
#colnames(nam_centromere) <- c("NAM_line", "Chr","Start","End","Size","100N gap number","13N gap number","Known gap size","Fully assembled (Y/n)","CentC","CRM2","CRM1","cinful-zeon","opie-ji","huck","prem1","grande", "All_TE","gene","Known sequences", "Other TE","Total gap","Unknown" )
colnames(nam_centromere) <- c("NAM_line", "Chr","Start","End","Size","100N gap number","13N gap number","Known gap size","Fully assembled (Y/n)","CentC","CRM","flip","prem1","cinful-zeon","gyma","doke","iteki","grande", "leviathan","centA", "uwum", "Copia", "all TE", "gene","Known sequences", "Other TE","Total gap","Unknown" )

head(nam_centromere)

nam_centromere_plot<-nam_centromere[ -c(23,25) ] # remove All_TE, Known sequences
nam_centromere_plot<-melt(nam_centromere_plot, id=c("NAM_line", "Chr","Start","End","Size","100N gap number","13N gap number","Known gap size","Fully assembled (Y/n)"))

head(nam_centromere_plot)
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
nam_centromere_plot$Chr <- factor(nam_centromere_plot$Chr,
                       levels = chr)
nam_centromere_plot$highlight[(nam_centromere_plot$`Fully assembled (Y/n)` == "Y")&(nam_centromere_plot$variable == "CentC")]  <- "*"

png("centromere_content_allchrs.png", units="in", width = 8, height = 7,res=300)
ggplot(nam_centromere_plot, aes(fill=variable, y=value, x=NAM_line)) + 
  geom_bar(position="stack", stat="identity",colour=NA) + 
  scale_y_continuous(labels=function(x)x/1000000)+ 
  ylab("Length of Components in active centromeres (Mb)") + 
  xlab("NAM lines") + 
  scale_fill_manual(values=c("#87a2ba", "#334863", "#eda64a","#666666","#ead09e","#604444","#2f4335","#1c8b82","#ab6661","#79538e","#a2a5b4","#dfc064","#90766b","#f4a3a4","#ffe2e2","#ec853f","#60965b")) +
  labs(fill = "Content") +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black")) 
dev.off()

pdf("centromere_content.pdf", width = 8, height = 9)
ggplot(nam_centromere_plot, aes(fill=variable, y=value, x=NAM_line)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_text(data = nam_centromere_plot,
            aes(label = highlight), 
            vjust = -1.2, size=9)+
  scale_y_continuous(labels=function(x)x/1000000)+ 
  ylab("Length of Components in active centromeres (Mb)") + 
  xlab("NAM lines") + 
  scale_fill_manual(values=c("#87a2ba", "#334863", "#eda64a","#666666","#ead09e","#604444","#2f4335","#1c8b82","#ab6661","#79538e","#a2a5b4","#dfc064","#90766b","#f4a3a4","#ffe2e2","#ec853f","#60965b")) +
  labs(fill = "Content") +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black")) +
  facet_grid(Chr~.) 
dev.off()
