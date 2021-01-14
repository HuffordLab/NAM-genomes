###NAM paper contribution--plotting
library("ggplot2")
library("stringr")
library("dplyr")
#install.packages("rgl", repos="http://R-Forge.R-project.org")
#BiocManager::install("karyoploteR")
library("karyoploteR")
library(data.table)

`%!in%` = Negate(`%in%`)

###ARRAY DATA. This is output from Rep_content_dat_prep.R
setwd("~/Desktop/NAM_contri/") 
rep_coords<- read.table("NAM_array_coords.tsv", header = T) ##All Array data from the Rep_content_dat_prep.R script
rep_coords<- rep_coords[rep_coords$V1 != "B73_AB10",] ##Removed Ab10
chr_num<- chr_num<- paste( "chr", 1:10, sep="")
comb_coords_chr_only<- rep_coords[rep_coords$Chr %in% chr_num,]

###REPEAT/TE data
path="~/Desktop/Dawe/TEST2"
all_gff<- dir(path, pattern="all_", recursive=TRUE, full.names=TRUE) ##all the summary files

lines1<- str_split(all_gff, pattern="_", simplify=T)[,2]
lines2<- str_split(lines1, pattern="[.]", simplify=T)[,1]
satel<- c("knob180", "TR-1", "CL569186.1", "CRM1", "repeat" , "(GGTTAG)n", "spacer", "CentC",lines2)
#remove the satellites-- this information should be in the array information already. smaller satellites 
# will be categorized as "other"

elem_list<- list()

for( i in 1:length(all_gff)){
  dat<-read.table(all_gff[i], sep="\t")
  new_dat_red<- dat[dat$V4 %!in% satel,]

  new_dat_red$V5<- NA

  for( k in 1:nrow(new_dat_red)){ ##removing nested hits
    bigger<- new_dat_red[new_dat_red$V1 == new_dat_red$V1[k] & new_dat_red$V2 <= new_dat_red$V2[k] & new_dat_red$V3 >=  new_dat_red$V3[k], ]
   if(nrow(bigger) > 1 ){new_dat_red$V5[k]<-"remove"}
  }

  new_dat_red2<- new_dat_red[is.na(new_dat_red$V5),] #grab the elements not tagged to be removed
  for(j in 2:nrow(new_dat_red2)){ ##removing overlaps
    if( new_dat_red2$V3[j-1] > new_dat_red2$V2[j]){new_dat_red2$V2[j] <- new_dat_red2$V3[j-1]+1} 
  }
    
  new_dat_red2$V5<- new_dat_red2$V3-new_dat_red2$V2 #this column should be empty, so overwrite with length
  new_dat_red3<-new_dat_red2[new_dat_red2$V5 >= 1,] #making sure they still have length after removing overlaps
  elem_list[[i]]<-new_dat_red3 #this list object has all TE's after edited and coordinate corrections for each line's knobs
  print(lines2[i])
}


#Find the biggest knobs
Biggest_knob<- as.data.frame(matrix(nrow=length(lines2)*10, ncol=ncol(rep_coords)))
big_knob_list<- list()
knob_type<-c("knob180", "TR-1")

for(i in 1:length(lines2)){
  lin<- lines2[i]
  elems<- elem_list[[i]]
  
  for(j in 1:length(chr_num)){
    lin_dat<- rep_coords[tolower(rep_coords$V1) == tolower(lin) & rep_coords$Chr == chr_num[j] & rep_coords$type %in% knob_type,]
    lin_dat_or<-  lin_dat[order(-lin_dat$Array.size..bp.),]
    Biggest_knob[(i-1)*10 + j,] <- lin_dat_or[1,] ##biggest knob

    tar_array_chr<-lin_dat_or$Chr[1]
    tar_array_start<-lin_dat_or$Start[1]
    tar_array_end<-lin_dat_or$End[1]
    
    #get arrays in that knob
    sub<- elems[elems$V1 == tar_array_chr & elems$V2 > tar_array_start & elems$V3 < tar_array_end, ]
    big_knob_list[[(i-1)*10 + j]] <- sub
  }
}


#Prep for plotting
#######
nam_knob_plot<- as.data.frame(matrix(nrow=0, ncol=5))
colnames(nam_knob_plot)<- c("NAM_line","array","variable","variable2", "value")

for(i in 1:length(big_knob_list)){
dat<-  big_knob_list[[i]]
not_tar<- c("Copia","Gypsy", "gene", satel)

dat$V6<- str_split(dat$V4, pattern="_", simplify=T)[,1]
dat$V7<- str_split(dat$V4, pattern="/", simplify=T)[,2]

TE_knobs<-  dat[dat$V7 %!in% not_tar, ]
TE_knobs_par=as.data.frame(matrix(nrow=1, ncol=3))
TE_knobs_par[1,1]<- "Other TE"
TE_knobs_par[1,2]<- "Other TE"
TE_knobs_par[1,3]<- sum(TE_knobs$V5)
colnames(TE_knobs_par)<- c("variable", "variable2", "value")

gene_knobs<-  dat[dat$V6 == "gene",]
gene_knobs%>% select(V6, V7, V5) %>%  group_by(V6,V7) %>% summarise_each(funs(sum))-> gene_knobs_par
colnames(gene_knobs_par)<- c("variable", "variable2","value")

Copia_knobs<-  dat[dat$V7 == "Copia",]
if(nrow(Copia_knobs > 0)){
  Copia_knobs$V6<- "Copia"
  Copia_knobs$V7<- NA
} 
#Copia_knobs_par=ddply(Copia_knobs,.(V6,V7),summarize, nLength=sum(V5))
Copia_knobs%>% select(V6, V7, V5) %>%  group_by(V6,V7) %>% summarise_each(funs(sum))->  Copia_knobs_par
colnames(Copia_knobs_par)<- c("variable", "variable2", "value")

Gypsy_knobs<-  dat[dat$V7 == "Gypsy",]
Gypsy_knobs %>% select(V6, V7, V5) %>%  group_by(V6,V7) %>% summarise_each(funs(sum))-> Gypsy_knobs_par
colnames(Gypsy_knobs_par)<- c("variable", "variable2","value")

TE_info<- rbind(TE_knobs_par, Copia_knobs_par, Gypsy_knobs_par, gene_knobs_par)

TE_info$array<- paste(Biggest_knob$V2[i], Biggest_knob$V3[i], sep="_")
TE_info$NAM_line<- Biggest_knob$V1[i]
TE_info_arr<- select(TE_info, c("NAM_line","array", "variable","variable2", "value"))

Rep_info<- as.data.frame(matrix(nrow=7, ncol=5))
colnames(Rep_info)<- c("NAM_line","array","variable", "variable2", "value")
Rep_info$NAM_line<- Biggest_knob$V1[i]
Rep_info$array<- paste(Biggest_knob$V2[i], Biggest_knob$V3[i], sep="_")

Rep_info$variable[1]<- "Total gap"
Rep_info$value[1]<- Biggest_knob$V6[i] * Biggest_knob$V10[i] 

Rep_info$variable[2]<- "rDNA"
Rep_info$value[2]<- Biggest_knob$V7[i]  * Biggest_knob$V11[i] 

Rep_info$variable[3]<- "CentC"
Rep_info$value[3]<- Biggest_knob$V7[i]  * Biggest_knob$V12[i] 

Rep_info$variable[4]<- "Subtelomeric Spacer"
Rep_info$value[4]<- Biggest_knob$V7[i]  * Biggest_knob$V13[i] 

Rep_info$variable[5]<- "Knob180"
Rep_info$value[5]<- Biggest_knob$V7[i]  * Biggest_knob$V14[i] 

Rep_info$variable[6]<- "TR-1"
Rep_info$value[6]<- Biggest_knob$V7[i]  * Biggest_knob$V15[i] 

Rep_info$variable[7]<- "Unknown"
Rep_info$value[7]<- Biggest_knob$V6[i]  - sum( Rep_info$value[1:6]) - sum(TE_info_arr$value)

total_info<- rbind(Rep_info, TE_info_arr)

nam_knob_plot<- rbind(nam_knob_plot, total_info)

}

nam_knob_plot$Chr<- str_split(nam_knob_plot$array, pattern="_", simplify=T)[,1]


nam_knob_plot %>% select(variable, variable2, value) %>% group_by(variable, variable2) %>% summarise_each(funs(sum)) -> check
#check=ddply(nam_knob_plot,.(variable, variable2) ,summarize, nLength=sum(value, na.rm = T))
check_or<- check[order(-check$value),]
rem<- c("CentC", "rDNA")
nam_knob_plot2<- nam_knob_plot[nam_knob_plot$variable %!in% rem,]

not_gypsy<-c("Knob180","TR-1", "Copia", "Other TE", "Unknown" , "Total gap" ,"Subtelomeric Spacer", "gene") #for variable2
keep_gypsy<- c("cinful", "flip","gyma","huck","xilon","doke") #stop at doke here. The total doke bp is 4328308.0, for variable
# and the next highest is prem1 with 2476141.0

nam_knob_plot2$unique_id<- paste(nam_knob_plot2$NAM_line, nam_knob_plot2$array, sep="_")

new_nam_knob_plot<- as.data.frame(matrix(nrow=0, ncol=ncol(nam_knob_plot2)))
for(i in 1:length(unique(nam_knob_plot2$unique_id))){
  sub<- nam_knob_plot2[nam_knob_plot2$unique_id == unique(nam_knob_plot2$unique_id)[i],]
  sub_notgyp<- sub[sub$variable2 %in% not_gypsy | is.na(sub$variable2) | sub$variable == "gene",]
  sub_gyp<- sub[sub$variable2 %!in% not_gypsy & !is.na( sub$variable2),]
  
  gyp_new<- as.data.frame(matrix(nrow=1, ncol=ncol(sub_gyp)))
  colnames(gyp_new)<- colnames(sub_gyp)
  sub_gyp2<- sub_gyp[sub_gyp$variable %in% keep_gypsy, ]
  sub_gyp3<- sub_gyp[sub_gyp$variable %!in% keep_gypsy, ]
  gyp_new$NAM_line[1]<- sub_gyp2$NAM_line[1]
  gyp_new$array[1]<- sub_gyp2$array[1]
  gyp_new$Chr[1]<- sub_gyp2$Chr[1]
  gyp_new$unique_id[1]<- sub_gyp2$unique_id[1]
  gyp_new$variable[1]<- "Other Gypsy"
  gyp_new$value[1]<- sum(sub_gyp3$value)
  
  mod_dat<-rbind(sub_notgyp, sub_gyp2,gyp_new)
    
  new_nam_knob_plot<- rbind(new_nam_knob_plot, mod_dat)
}


Biggest_knob$unique_id<- paste(Biggest_knob$V1, Biggest_knob$V2, Biggest_knob$V3, sep="_")
new_nam_knob_plot_highlight<- as.data.frame(matrix(nrow=0, ncol=ncol(new_nam_knob_plot)+1))

for(i in 1:length(unique(Biggest_knob$unique_id))){
  tar<- unique(Biggest_knob$unique_id)[i]
  
  hili<- Biggest_knob[Biggest_knob$unique_id %in% tar,]$V9
  sub<- new_nam_knob_plot[new_nam_knob_plot$unique_id %in% tar ,]
  sub$highlight<-NA
  if( !is.na(hili)){
  sub[sub$variable=="Unknown",]$highlight<- "*"
  sub[sub$variable!="Unknown",]$highlight<- " "
  } else{sub$highlight<- " " }
  new_nam_knob_plot_highlight<- rbind(new_nam_knob_plot_highlight, sub)
  
}

new_nam_knob_plot_highlight$variable[new_nam_knob_plot_highlight$variable == "Subtelomeric Spacer"]<- "subtelomeric repeat"
new_nam_knob_plot_highlight$variable[new_nam_knob_plot_highlight$variable == "cinful"]<- "cinful-zeon"
new_nam_knob_plot_highlight$variable[new_nam_knob_plot_highlight$variable == "Unknown"]<- "un-annotated"
new_nam_knob_plot_highlight$variable[new_nam_knob_plot_highlight$variable == "Knob180"]<- "knob 180"
new_nam_knob_plot_highlight$variable[new_nam_knob_plot_highlight$variable == "Other Gypsy"]<- "other gypsy"
new_nam_knob_plot_highlight$variable[new_nam_knob_plot_highlight$variable == "Other TE"]<- "other TE"
new_nam_knob_plot_highlight$variable[new_nam_knob_plot_highlight$variable == "Total gap"]<- "gaps"
new_nam_knob_plot_highlight$variable[new_nam_knob_plot_highlight$variable == "Copia"]<- "copia"


new_nam_knob_plot_highlight$variable<-as.factor(new_nam_knob_plot_highlight$variable)
new_nam_knob_plot_highlight$variable= factor(new_nam_knob_plot_highlight$variable, levels(new_nam_knob_plot_highlight$variable)[c(9, 13, 12, 1, 4, 7, 8, 15, 3, 10, 2, 6, 11, 5, 14)])

new_nam_knob_plot_highlight$Chr<-as.factor(new_nam_knob_plot_highlight$Chr)
new_nam_knob_plot_highlight$Chr= factor(new_nam_knob_plot_highlight$Chr, levels(new_nam_knob_plot_highlight$Chr)[c(1, 3:10, 2)])

new_nam_knob_plot_highlight$NAM_line<-as.factor(new_nam_knob_plot_highlight$NAM_line)
new_nam_knob_plot_highlight$NAM_line= factor(new_nam_knob_plot_highlight$NAM_line, levels(new_nam_knob_plot_highlight$NAM_line)[c(1:2, 15:16, 19, 22:23, 17:18, 25, 11, 24, 12, 9:10, 3:8, 14, 13,20:21,26)])

new_nam_knob_plot_highlight<- new_nam_knob_plot_highlight[!is.na(new_nam_knob_plot_highlight$Chr),]

x_cols <- c(rep(c("goldenrod1"),1),rep(c("gray47"),3),rep(c("royalblue"),6),rep(c("orchid"),1),rep(c("orangered"),2),rep(c("limegreen"),13))

knob_plot<-ggplot(new_nam_knob_plot_highlight, aes(fill=variable, y=value, x=NAM_line)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_text(data = new_nam_knob_plot_highlight,
           aes(label = highlight),vjust = -1.2,
             size=9)+
  scale_y_continuous(labels=function(x)x/1000000)+ 
  ylab("Length of Components in Largest Knob Arrays (Mb)") + 
  scale_fill_manual(values=colors) +
  labs(fill = "Content") +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=11,angle = 40, hjust = 1, colour= x_cols),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=15,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +facet_grid(Chr~., scales="free_y") 
knob_plot

colors<-c( "#1c8b82",  "#ab6661", "#eda64a","#666666","#ead09e","#604444","#2f4335",
           "#79538e","#a2a5b4","#dfc064","#90766b","#f4a3a4","#ffe2e2","#ec853f","#60965b")

ggsave( plot=knob_plot, "biggestknobs_resize_v2.png", dpi=200, width=8, height=9)

write.table(new_nam_knob_plot_highlight, "new_nam_knob_plot_highlight.txt", quote=F, row.names = F, sep="\t")
####

#part B
#setwd("~/Desktop/NAM_contri/") 
new_nam_knob_plot_highlight2<- read.table(file="new_nam_knob_plot_highlight.txt", header = T, sep="\t")
top5<- levels(new_nam_knob_plot_highlight$variable)[1:7]
top5<- c( top5[1:2], top5[4:7])

high_array_dat<- new_nam_knob_plot_highlight[new_nam_knob_plot_highlight$NAM_line == "B73" & new_nam_knob_plot_highlight$Chr == "chr5",]
Biggest_knob_dat<- Biggest_knob[Biggest_knob$V1 == "B73" & Biggest_knob$V2 == "chr5",]
#using big knob on chr 1 as example. That will be 1 in the big_knob_list data
big_knob_elems<- big_knob_list[[5]]
big_knob_elems$V6<- str_split(big_knob_elems$V4, pattern="_", simplify=T)[,1]
big_knob_elems$V7<- str_split(big_knob_elems$V4, pattern="/", simplify=T)[,2]


pp <- getDefaultPlotParams(plot.type = 2)
pp$data2height <- 90
pp$ideogramheight <- 0.5
pp$data1inmargin<-0

genome <- data.frame("line" = c('chr5'), "Start" = c(199000000), "End" = c(204000000))
Knob180<-big_knob_elems[big_knob_elems$V6 == "knob180", ]
TR<-big_knob_elems[big_knob_elems$V6 == "TR-1", ]
cinful<-big_knob_elems[big_knob_elems$V6 == "cinful", ]
flip<-big_knob_elems[big_knob_elems$V6 == "flip", ]
gyma<-big_knob_elems[big_knob_elems$V6 == "gyma", ]
#huck<-big_knob_elems[big_knob_elems$V6 == "huck", ]

dev.off()
dev.new()
chromosome <- plotKaryotype(genome = genome, plot.type=2, plot.params = pp)
kpRect(chromosome, data=toGRanges(Knob180),x0=Knob180$V2,x1=Knob180$V3,y0=0, y1=1,lwd=0.00,r0=0.60, r1=0.70,col="#1c8b82",border=NA)
kpRect(chromosome, data=toGRanges(TR),x0=TR$V2,x1=TR$V3,y0=0, y1=1,lwd=0.00,r0=0.48, r1=0.58,col="#ab6661",border=NA,clipping = TRUE)
kpRect(chromosome, data=toGRanges(cinful),x0=cinful$V2,x1=cinful$V3,y0=0, y1=1,lwd=0.00,r0=0.36, r1=0.46,col='#666666',border=NA,clipping = TRUE)
kpRect(chromosome, data=toGRanges(flip),x0=flip$V2,x1=flip$V3,y0=0, y1=1,lwd=0.00,r0=0.24, r1=0.34,col="#ead09e",border=NA,clipping = TRUE)
kpRect(chromosome, data=toGRanges(gyma),x0=gyma$V2,x1=gyma$V3,y0=0, y1=1,lwd=0.00,r0=0.12, r1=0.22,col="#604444",border=NA,clipping = TRUE)
#kpRect(chromosome, data=toGRanges(huck),x0=huck$V2,x1=huck$V3,y0=0, y1=1,lwd=0.00,r0=0, r1=0.1,col="#2f4335",border=NA,clipping = TRUE)
kpAddBaseNumbers(chromosome, tick.dist = 1000000, add.units = FALSE,digits=1,cex=0,tick.len = 15,minor.tick.dist = 1000000, minor.tick.len = 5,clipping=TRUE)
png("B73_4_chr5.png", width = 7, height = 2, units = 'in', res = 500)
