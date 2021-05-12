library(tidyverse)
library(ggplot2)

#Bring in expression data
##Read in the raw TPM for each genome and the key to the gene ids
GeneIDKey<-read_delim("~/Desktop/Grad School/NAM_FloweringTime_Project/Expression_data/allgeneids_wB73v5.txt", delim = " ")
colnames(GeneIDKey)<-c("V3","B97","CML103","CML228","CML247",
                                   "CML277","CML322","CML333","CML52",
                                   "CML69","HP301","Il14H","Ki11","Ki3",
                                   "Ky21","M162W","M37W","Mo18W","Ms71",
                                   "NC350","NC358","Oh43","Oh7B","P39",
                                   "Tx303","Tzi8","B73")

#Make a master TPM dataframe
AllTPM<-read_tsv("~/Desktop/Grad School/NAM_FloweringTime_Project/Expression_data/allTPM_splitdescriptions.txt", col_names = TRUE)
AllTPM<-mutate(AllTPM, gene=substring(AllTPM$samples,9,22))
GeneIDKey %>% separate(B97, c("B97_1","B97_2","B97_3"), sep = ":")

NAMnames<-c("B97","CML103","CML228","CML247",
            "CML277","CML322","CML333","CML52",
            "CML69","HP301","Il14H","Ki11","Ki3",
            "Ky21","M162W","M37W","Mo18W","Ms71",
            "NC350","NC358","Oh43","Oh7B","P39",
            "Tx303","Tzi8","B73")
#This splits up the Geneid names with the colons
wide_geneIDKey<-select(GeneIDKey, V3)
for(i in 1:length(NAMnames)){
  temp<-GeneIDKey %>% separate(NAMnames[i], c(paste(NAMnames[i],1,sep = "_"),
                                             paste(NAMnames[i],2,sep = "_"),
                                             paste(NAMnames[i],3,sep = "_")), sep = ":") %>% select(contains("_"))
  wide_geneIDKey<-cbind(wide_geneIDKey, temp)
}
#these two steps get the V3 ID given a gene from TPM DF
temp<-colnames(wide_geneIDKey[grepl(AllTPM$gene[1],wide_geneIDKey)])
filter(wide_geneIDKey, get(temp)==AllTPM$gene[1]) %>% select(V3)
AllTPM$V3<-NA #Creates a new column 
for(i in 1:nrow(AllTPM)){
  temp<-colnames(wide_geneIDKey[grepl(AllTPM$gene[i],wide_geneIDKey)])
  AllTPM$V3[i]<-filter(wide_geneIDKey, get(temp)==AllTPM$gene[i]) %>% select(V3) %>% pull()
}

#Add TPM to the haplotype dataframes (requires the FT_Indel_Analysis.Rmd to be run)
#Test it out
#V5 IDs are in wide_geneIDKey$B73_1
#Putting the significant findings into DF 
sig_indels<-tibble(gene=c("Zm00001e009095","Zm00001e035953","Zm00001e001462","Zm00001e009095","Zm00001e022997" ),
                   indel=c("ins23110","del76968","del3270","ins23110","ins45095" ))


for(i in 1:nrow(sig_indels)){
  ls(pattern = sig_indels$gene[i]) %>% print()
}
sig_indels<-add_column(sig_indels, df=c("G.Zm00001e009095 chr2:145013763-145019458","G.Zm00001e035953 chr7:183947136-183980571",
                                        "P.Zm00001e001462 chr1:49874953-49879953","P.Zm00001e009095 chr2:145008763-145013763",
                                        "P.Zm00001e022997 chr4:145822771-145827771"))


#trial<-ls(pattern = wide_geneIDKey$B73_1[1]) #get the names for the haplotype DFs
#get(trial[1])
#wide_geneIDKey$V3[1] #Gives us the shared V3 ID across genomes
#trial1<-inner_join(filter(AllTPM, V3 == wide_geneIDKey$V3[1]),
#                   (grep(wide_geneIDKey$B73_1[1], sig_indels$df, value = TRUE) %>% get()), 
#                   by="Genome")
#trial1
#to look at relationship between flowering time and TPM
#ggplot(data = trial1) +
#  geom_point(aes(x=Flower_Time, y=TPM))
#to look at relationship between indels and TPM
#ggplot(data = trial1) +
#  geom_jitter(aes(x=as.factor(ins117310),y=TPM),width = 0.2)+
#  geom_boxplot(aes(x=as.factor(ins117310), y=TPM),color="blue", alpha=0.5) 
  #geom_boxplot(aes(x=as.factor(INS34657), y=TPM),color="orange", alpha=0.5)

#inner_join(filter(AllTPM, V3 == wide_geneIDKey$V3[1]),
#           (grep(wide_geneIDKey$B73_1[1], sig_indels$df, value = TRUE) %>% get()), 
#           by="Genome") %>%
#  ggplot()+
#  geom_jitter(aes(x=as.factor(ins117310),y=TPM),width = 0.2)+
#  geom_boxplot(aes(x=as.factor(ins117310), y=TPM),color="blue", alpha=0.5) 

#make_haplo_plot<-function(g){
#  t<-inner_join(filter(AllTPM, V3 == wide_geneIDKey$V3[g]),
#                (grep(wide_geneIDKey$B73_1[g], sig_indels$df, value = TRUE) %>% get()),
#                by="Genome")
#  i<-which(sig_indels$gene %in% c(wide_geneIDKey$B73_1[g]))
#  for(j in 1:length(i)){
#    indel<-sig_indels$indel[i[j]]
#    plt<-ggplot(data=t)+
#      geom_jitter(aes_string(x=indel,y="TPM",group=indel),width = 0.2)+
#      geom_boxplot(aes_string(x=indel,y="TPM",group=indel), color="blue",alpha=0.5)+
#      ggtitle(paste(wide_geneIDKey$B73_1[g],indel,t$Region_Type[1],sep=" "))
#    return(plt)
# }
#}

#make_haplo_plot(1)
#make_haplo_plot(2)#still doesn't work if there are multiple sig indels

#Significant indels
#add V3 id names as separate column to help
sig_indels<-add_column(sig_indels,V3=c("HQ003893.1","Zm00001d022590_T004","GRMZM2G124532_T01","HQ003893.1","GRMZM2G020081_T01"))
#Trial run of making the figures

#inner join of the TPM values with the V3 id name with the dataframe from sig indels, joined by genome
t<-inner_join(filter(AllTPM, V3 == sig_indels$V3[1]),(get(sig_indels$df[1])),by="Genome")
#indel gets us the significant indel for that row
indel<-sig_indels$indel[1]
indel
#makes the box plot
#plt<-
ggplot(data = t, aes(x=Tissue, y=TPM, fill=as.factor(ins23110)))+ 
  geom_boxplot()

test<-t
test$ins23110<-test$ins23110 %>% as.factor()
ggplot(data = test, aes(x=Tissue))+
  geom_jitter(aes_string(y="TPM", shape=indel),width = 0.2, alpha=0.5)+
  geom_boxplot(aes_string(y="TPM",fill=indel),alpha=0.75)+
  facet_wrap(~ins23110)+
  theme(axis.text.x=element_text(angle=60))

test$Tissue<-factor(test$Tissue, levels=c("V11_base","V11_middle","V11_tip","V18_ear","V18_tassel","R1_anther"),ordered=TRUE)
ggplot(data = test, aes_string(x=indel))+
  geom_jitter(aes_string(y="TPM", color=indel),width = 0.2, alpha=0.5)+
  geom_boxplot(aes_string(y="TPM"), alpha=0.75)+
  facet_wrap(~Tissue, scales = "free")

test[,indel]

#makes the flowering time plot
ggplot(data=test)+
  geom_point(aes_string(x="Flower_Time",y="TPM",color=indel))+
  geom_smooth(aes_string(x="Flower_Time",y="TPM",group = indel, linetype=indel),color="black")+
  facet_wrap(~Tissue)+
  ggtitle(paste(sig_indels$gene[5],"expression v DTA",sep = " "))

#Loop to make the figures

#t<-inner_join(filter(AllTPM, V3 == sig_indels$V3[2]),(get(sig_indels$df[2])),by="Genome")
for(i in 1:nrow(sig_indels)){
  #inner join of the TPM values with the V3 id name with the dataframe from sig indels, joined by genome
  t<-inner_join(filter(AllTPM, V3 == sig_indels$V3[i]),(get(sig_indels$df[i])),by="Genome")
  #indel gets us the significant indel for that row
  indel<-sig_indels$indel[i]
  indel
  t[,indel]<- t[,indel] %>% pull() %>% as.factor()
  t$Tissue<-factor(t$Tissue, levels=c("V11_base","V11_middle","V11_tip","V18_ear","V18_tassel","R1_anther"),ordered=TRUE)
  #makes the box plot
  plt<-ggplot(data = t,aes_string(x=indel))+
    geom_jitter(aes_string(y="TPM", color=indel),width = 0.2, alpha=0.5)+
    geom_boxplot(aes_string(y="TPM"), alpha=0.75)+
    facet_wrap(~Tissue, scales = "free")
    ggtitle(paste(sig_indels$gene[i],indel,t$Region_Type[1],sep=" "))
  ggsave(plt,filename = paste("Expression_data/Plots/",sig_indels$gene[i],indel,t$Region_Type[1],".png",sep=""))
  #makes the flowering time plot
  #plt<-ggplot(data=t)+
   # geom_point(aes_string(x="Flower_Time",y="TPM",color=indel, group=indel))+
    #ggtitle(paste(sig_indels$gene[i],"TPM v DTA",t$Region_Type[1],sep = " "))
  #ggsave(plt,filename = paste("Expression_data/Plots/",sig_indels$gene[i],indel,"TPM v DTA",t$Region_Type[1],".png",sep = ""))
}

###Zm00001d022590_T004 doesn't have a P39 gene id and the sig indel for it is only
###in P39. Check if this is an issue with the blast or join.
###
### Did get a blast match in P39
#[snodgras@condo2017 blast_results]$ grep Zm00001d022590_T004 tophit_out_P39_nostrd.bed 
#chr7	177903103	177903894	Zm00001d022590_T004
#Does not have a gene id in the P39 coordinate system, was not annotated in the GFF

###Making box plots for all the indels tested that were not significant
#for each dataframe, 
  #pull out the colnames of the SVs
  #add to a dataframe: df = dataframe name, indel= SV name, gene=V5name, V3=V3ID

#`G.no_evd_annotation chr1:275855389-275861912`[,5:ncol(`G.no_evd_annotation chr1:275855389-275861912`)] %>% colnames()

placeholder<-ls(pattern = ":")
#SV<-get(placeholder[1]) %>% select(-c(Region,Region_Type,Genome,Flower_Time)) %>% colnames()
#for(j in 1:length(SV)){print(SV[j])}

df<-c()
indel<-c()
for(i in 1:length(placeholder)){
  SV<-get(placeholder[i]) %>% select(-c(Region,Region_Type,Genome,Flower_Time)) %>% colnames()
  for(j in 1:length(SV)){
    df<-c(df,placeholder[i])
    indel<-c(indel,SV[j])
  }
}
library(stringi)
gene<-c()
for(k in 1:length(df)){
  temp<-stri_split_fixed(df[k], pattern=" ")[[1]][1]
  gene<-c(gene, stri_split_fixed(temp,pattern = ".")[[1]][2])
}

non_sig_indels<-tibble(df=df,indel=indel,gene=gene)

#Need to add in V3 id
V3<-c("GRMZM2G157727_T01","GRMZM2G181028_T01","GRMZM2G124532_T01","GRMZM2G092174_T01","GRMZM2G057935_T01","GRMZM2G129889_T01","GRMZM5G861678_T01","GRMZM2G107101_T01","GRMZM5G844173_T01","GRMZM2G067702_T01","GRMZM2G014902_T01","GRMZM2G107945_T01","GRMZM2G106363_T01","EF114229.2","GRMZM2G474769_T01","Zm00001d022590_T004" ,"HQ003893.1","EU952116.1","GRMZM2G020081_T01","GRMZM2G405368_T01","GRMZM2G381691_T01","GRMZM2G179264_T01","GRMZM2G011357_T01","AF166527.1","NM_001254783.1","GRMZM2G144744_T01","GRMZM2G024973_T01","GRMZM2G078798_T01","GRMZM2G017087_T01","GRMZM2G067921_T01","GRMZM2G148693_T01" ,"GRMZM2G072582_T01" ,"GRMZM2G156079_T01","GRMZM2G098813_T01","GRMZM2G180190_T01","GRMZM2G032339_T01","GRMZM2G171365_T01","GRMZM2G700665_T01","GRMZM2G004483","GRMZM2G021777")
V5<-c("no_evd_annotation","Zm00001e005499","Zm00001e001462","Zm00001e038324","Zm00001e005702","Zm00001e012030","Zm00001e028329","Zm00001e025967","Zm00001e016671","Zm00001e028258","Zm00001e021724","Zm00001e011193","Zm00001e024573","Zm00001e032705","Zm00001e021724","Zm00001e035953","Zm00001e009095","Zm00001e038444","Zm00001e022997","Zm00001e036812","Zm00001e040482","Zm00001e027418","Zm00001e004751","Zm00001e018725","Zm00001e037501","Zm00001e005412","Zm00001e012213","Zm00001e031130","Zm00001e005545","Zm00001e035978",'Zm00001e032703',"Zm00001e032703","Zm00001e022651","Zm00001e041565","ZFL2","Zm00001e005705","Zm00001e039048","Zm00001e027593","ZmCCT9","ZmCOL3")
V3toV5<-tibble(V3=V3,gene=V5)

non_sig_indels<-inner_join(non_sig_indels,V3toV5)
#These commands removes joined rows that have wrong names because some V3 ids have identical V5 ids
non_sig_indels<-filter(non_sig_indels, !V3 == "GRMZM2G148693_T01" & !df == "G.Zm00001e032703 chr7:2156721-2168013" )
non_sig_indels<-filter(non_sig_indels, !V3 == "GRMZM2G072582_T01" & !df == "G.Zm00001e032703 chr2:242277539-242286715" )

non_sig_indels<-filter(non_sig_indels, !indel %in% sig_indels$indel)

for(i in 1:nrow(non_sig_indels)){
  #inner join of the TPM values with the V3 id name with the dataframe from sig indels, joined by genome
  t<-inner_join(filter(AllTPM, V3 == non_sig_indels$V3[i]),(get(non_sig_indels$df[i])),by="Genome")
  #indel gets us the significant indel for that row
  indel<-non_sig_indels$indel[i]
  indel
  t[,indel]<- t[,indel] %>% pull() %>% as.factor()
  t$Tissue<-factor(t$Tissue, levels=c("V11_base","V11_middle","V11_tip","V18_ear","V18_tassel","R1_anther"),ordered=TRUE)
  #makes the box plot
  plt<-ggplot(data = t,aes_string(x=indel))+
    geom_jitter(aes_string(y="TPM", color=indel),width = 0.2, alpha=0.5)+
    geom_boxplot(aes_string(y="TPM"), alpha=0.75)+
    facet_wrap(~Tissue, scales = "free")
  ggtitle(paste(non_sig_indels$gene[i],indel,t$Region_Type[1],sep=" "))
  ggsave(plt,filename = paste("Expression_data/Plots/",non_sig_indels$gene[i],indel,t$Region_Type[1],".png",sep=""))
}

###########################################################
###########################################################
#t-tests for unequal sample sizes

#g<-which(wide_geneIDKey$B73_1 %in% c("Zm00001e039048")) #get index position in the widegene key
#wide_geneIDKey[g,c(1,77)] #verify
#t<-inner_join(filter(AllTPM, V3 == wide_geneIDKey$V3[g]),
#              (grep(wide_geneIDKey$B73_1[g], sig_indels$df, value = TRUE) %>% get()),
#              by="Genome")
t<-inner_join(filter(AllTPM, V3 == sig_indels$V3[1]),(get(sig_indels$df[1])),by="Genome")

t<-t %>% select(c(Tissue, TPM, sig_indels$indel[1])) #subsets important groups

CG<-filter(t, t[,3] == 0) #sets those without indel to control group
TG<-filter(t,t[,3] == 1) #sets those with indel to treatment group

tissues<-t %>% select(Tissue) %>% unique() %>% pull()

CG_tis<-filter(CG, Tissue == tissues[1])
TG_tis<-filter(TG, Tissue == tissues[1])
#var.test(CG$TPM,TG$TPM)#checks to see if the variance is equal
var.test(CG_tis$TPM,TG_tis$TPM)$p.value
agep<-t.test(CG_tis$TPM,TG_tis$TPM, var.equal = TRUE)
agep$estimate
agep$p.value

#Making a loop for significant Indels
tissues<-t %>% select(Tissue) %>% unique() %>% pull()
est.x<-c()
est.y<-c()
p<-c()
tis<-c()
ind<-c()
for(i in 1:nrow(sig_indels)){
  t<-inner_join(filter(AllTPM, V3 == sig_indels$V3[i]),(get(sig_indels$df[i])),by="Genome")
  t<-t %>% select(c(Tissue, TPM, sig_indels$indel[i]))
  for(j in 1:length(tissues)){
    CG<-filter(t, t[,3] == 0 & Tissue == tissues[j])
    TG<-filter(t,t[,3] == 1 & Tissue == tissues[j])
    if(nrow(CG) > 1 & nrow(TG) > 1 & sum(CG$TPM,TG$TPM) != 0){
      if(var.test(CG$TPM,TG$TPM)$p.value <= 0.05){
        agep<-t.test(CG$TPM, TG$TPM, var.equal = FALSE)
        est.x<-c(est.x, agep$estimate[1])
        est.y<-c(est.y, agep$estimate[2])
        p<-c(p,agep$p.value)
        tis<-c(tis, tissues[j])
        ind<-c(ind, sig_indels$indel[i])
      }
      else{
        agep<-t.test(CG$TPM, TG$TPM, var.equal = TRUE)
        est.x<-c(est.x, agep$estimate[1])
        est.y<-c(est.y, agep$estimate[2])
        p<-c(p,agep$p.value)
        tis<-c(tis, tissues[j])
        ind<-c(ind, sig_indels$indel[i])
      }
    }
    else{
      est.x<-c(est.x, NA)
      est.y<-c(est.y, NA)
      p<-c(p,NA)
      tis<-c(tis, tissues[j])
      ind<-c(ind, sig_indels$indel[i]) 
    }
  }
}
#Loop to test for non-sigificant indels
for(i in 1:nrow(non_sig_indels)){
  t<-inner_join(filter(AllTPM, V3 == non_sig_indels$V3[i]),(get(non_sig_indels$df[i])),by="Genome")
  t<-t %>% select(c(Tissue, TPM, non_sig_indels$indel[i]))
  for(j in 1:length(tissues)){
    CG<-filter(t, t[,3] == 0 & Tissue == tissues[j])
    TG<-filter(t,t[,3] == 1 & Tissue == tissues[j])
    if(nrow(CG) > 1 & nrow(TG) > 1 & sum(CG$TPM,TG$TPM) != 0 ){
      if(var.test(CG$TPM,TG$TPM)$p.value <= 0.05){
        agep<-t.test(CG$TPM, TG$TPM, var.equal = FALSE)
        est.x<-c(est.x, agep$estimate[1])
        est.y<-c(est.y, agep$estimate[2])
        p<-c(p,agep$p.value)
        tis<-c(tis, tissues[j])
        ind<-c(ind, non_sig_indels$indel[i])
      }
      else{
        agep<-t.test(CG$TPM, TG$TPM, var.equal = TRUE)
        est.x<-c(est.x, agep$estimate[1])
        est.y<-c(est.y, agep$estimate[2])
        p<-c(p,agep$p.value)
        tis<-c(tis, tissues[j])
        ind<-c(ind, non_sig_indels$indel[i])
      }
    }
    else{
      est.x<-c(est.x, NA)
      est.y<-c(est.y, NA)
      p<-c(p,NA)
      tis<-c(tis, tissues[j])
      ind<-c(ind, non_sig_indels$indel[i]) 
    }
  }
}

#create a table with all the test results
test<-tibble(CG_est=est.x,TG_est=est.y,p_val=p,Tissue=tis,indel=ind)
0.05/372 #Bonferroni cut off
#Which of the significant indels meet the Bonferroni cut off?
filter(test, p_val <= 0.0001344086 & indel %in% sig_indels$indel)
#Which of the insignificant indels meet the Bonferroni cut off
filter(test, p_val <= 0.0001344086 & indel %in% non_sig_indels$indel) %>% print(n=100)

###Calculate Benjamini-Hochberg values for test values
test<-arrange(test, p_val) #arrange p-values smallest to largest
#Find the largest integer k so that p(k) <= k*alpha/m
BH_05<-c()
for(k in 1:nrow(test)){
  temp<-(k*0.05)/nrow(test)
  BH_05<-c(BH_05,temp)
}
test<-add_column(test,BH_05)
test
temp<-filter(test, p_val <= test$BH_05) %>% print(n=30)

#Find the number of lines with a given indel
filter(non_sig_indels, indel == "del14307") %>% select(df) %>% pull() %>% get() %>% filter(del14307 == 1) %>% nrow()

n<-c()
for(i in 1:nrow(temp)){
  ind<-temp$indel[i]
  if(TRUE %in% grepl(ind, non_sig_indels)){
    j<-filter(non_sig_indels, indel == ind) %>% select(df)%>% pull() %>% get()%>%select(c(Genome,ind))
    n<-c(n,filter(j, j[,2] == 1) %>% nrow())
  }
  else{
    j<-filter(sig_indels, indel == ind) %>% select(df)%>% pull() %>% get()%>%select(c(Genome,ind))
    n<-c(n,filter(j, j[,2] == 1) %>% nrow())
  }
  remove(j)
}
temp<-add_column(temp, Genome_Number = n)
#Genome number is the number of genomes that have this indel
filter(temp, Genome_Number > 1)

#Try log2 difference
#Info came from http://rstudio-pubs-static.s3.amazonaws.com/13988_bb11d85b79b2436280de434988558140.html
temp<-mutate(temp, fold_change= TG_est/CG_est,
       log2fc = log2(fold_change))

#Trying to figure out which genomes are the extremes for potential allelic series
t<-inner_join(filter(AllTPM, V3 == "NM_001254783.1"),`P.Zm00001e037501 chr9:101039144-101044144`,by="Genome")
#indel gets us the significant indel for that row
indel<-"ins88899" #"ins88900"
#indel
t[,indel]<- t[,indel] %>% pull() %>% as.factor()
t$Tissue<-factor(t$Tissue, levels=c("V11_base","V11_middle","V11_tip","V18_ear","V18_tassel","R1_anther"),ordered=TRUE)
#makes the box plot
ggplot(data = t,aes_string(x=indel))+
  geom_jitter(aes_string(y="TPM", color=indel),width = 0.2, alpha=0.5)+
  geom_boxplot(aes_string(y="TPM"), alpha=0.75)+
  geom_text(aes(y=TPM,label=Genome),size=3,alpha=0.6,position = "jitter")+
  facet_wrap(~Tissue, scales = "free")+
  ggtitle(paste("Zm00001e037501",indel,t$Region_Type[1],sep=" "))

temp<-add_column(temp, Gene=c("Zm00001e005545","Zm00001e012030","Zm00001e012030","Zm00001e027593",
                              "Zm00001e001462","Zm00001e035978","Zm00001e012030","no_evd_annotation",
                              "Zm00001e028258","no_evd_annotation","Zm00001e035953","no_evd_annotation",
                              "Zm00001e037501","Zm00001e038444","Zm00001e005702","Zm00001e022997",
                              "Zm00001e009095","no_evd_annotation","Zm00001e001462","Zm00001e009095",
                              "Zm00001e022997","Zm00001e022997","Zm00001e022997","Zm00001e037501"))
write_csv(temp,"Expression_data/sig_expression_test_results.csv")

#Try doing the haplotype thing
GL15_promoter<-select(`P.Zm00001e037501 chr9:101039144-101044144`, c(Region,Region_Type,Genome,Flower_Time))
GL15_promoter<-add_column(GL15_promoter, Indel_haplotype = c("900","900","900","none","none","899","none","900","none","none","both","both","both","none","none","899","none","none","none","none","none","900","both","none","both"))
#### Zm00028a052023 This might be the IL14H id for GL15
t<-inner_join(filter(AllTPM, V3 == "NM_001254783.1"),GL15_promoter,by="Genome")
t[,"Indel_haplotype"]<-factor(t$Indel_haplotype, levels = c("none","899","900","both"))
t$Tissue<-factor(t$Tissue, levels=c("V11_base","V11_middle","V11_tip","V18_ear","V18_tassel","R1_anther"),ordered=TRUE)
t<-add_row(t, Length = c(3600,3600,3600,3600,3600,3600,3600,3600,3600,3600,3600,3600), 
        samples=c("ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600",
                  "ID=gene:Zm00028a052023;biotype=protein_coding;logic_name=mikado_gene;length=3600"),
        Tissue_ID=c("cnts_Il14H_R1_anther_MN13081","cnts_Il14H_R1_anther_MN13082","cnts_Il14H_V11_base_MN13031",
                    "cnts_Il14H_V11_base_MN13032","cnts_Il14H_V11_middle_MN13041","cnts_Il14H_V11_middle_MN13042",
                    "cnts_Il14H_V11_tip_MN13051","cnts_Il14H_V11_tip_MN13052","cnts_Il14H_V18_ear_MN13071",
                    "cnts_Il14H_V18_ear_MN13072","cnts_Il14H_V18_tassel_MN13061","cnts_Il14H_V18_tassel_MN13062"),
        TPM=c(0.0342499404944243,0.044291267616968,0.462043383239728,0.220325494944133,5.33857386992319,10.1484701264353,
              5.60204186614248,4.61706409331874,1.0180143101381,0.800205253092941,0.596681443059441,0.468410744572754),
        Genome=c("IL14H","IL14H","IL14H","IL14H","IL14H","IL14H","IL14H","IL14H","IL14H","IL14H","IL14H","IL14H"),
        Tissue=c("R1_anther","R1_anther","V11_base","V11_base","V11_middle","V11_middle","V11_tip","V11_tip",
                 "V18_ear","V18_ear","V18_tassel","V18_tassel"),
        Replicate=c("MN13081","MN13082","MN13031","MN13032","MN13041","MN13042","MN13051","MN13052",
                    "MN13071","MN13072","MN13061","MN13062"),
        gene=c("Zm00028a052023","Zm00028a052023","Zm00028a052023","Zm00028a052023","Zm00028a052023",
               "Zm00028a052023","Zm00028a052023","Zm00028a052023","Zm00028a052023","Zm00028a052023",
               "Zm00028a052023","Zm00028a052023"),
        V3=c("NM_001254783.1","NM_001254783.1","NM_001254783.1","NM_001254783.1","NM_001254783.1","NM_001254783.1",
             "NM_001254783.1","NM_001254783.1","NM_001254783.1","NM_001254783.1","NM_001254783.1","NM_001254783.1"),
        Region=c("Zm00001e037501 chr9:101039144-101044144","Zm00001e037501 chr9:101039144-101044144",
                 "Zm00001e037501 chr9:101039144-101044144","Zm00001e037501 chr9:101039144-101044144",
                 "Zm00001e037501 chr9:101039144-101044144","Zm00001e037501 chr9:101039144-101044144",
                 "Zm00001e037501 chr9:101039144-101044144","Zm00001e037501 chr9:101039144-101044144",
                 "Zm00001e037501 chr9:101039144-101044144","Zm00001e037501 chr9:101039144-101044144",
                 "Zm00001e037501 chr9:101039144-101044144","Zm00001e037501 chr9:101039144-101044144"),
        Region_Type=c("Promoter","Promoter","Promoter","Promoter","Promoter","Promoter","Promoter",
                      "Promoter","Promoter","Promoter","Promoter","Promoter"),
        Indel_haplotype=c("900","900","900","900","900","900","900","900","900","900","900","900"))
#makes the box plot
ggplot(data = t,aes(x=Indel_haplotype))+
  geom_jitter(aes(y=TPM, color=Indel_haplotype),width = 0.2, alpha=0.5)+
  geom_boxplot(aes(y=TPM), alpha=0.75)+
  geom_text(aes(y=TPM,label=Genome),size=3,alpha=0.6,position = "jitter")+
  facet_wrap(~Tissue, scales = "free")+
  ggtitle(paste("Zm00001e037501",t$Region_Type[1],sep=" "))

t %>% filter(Tissue %in% c("V11_base","V11_middle","V11_tip"))%>%ggplot(aes(x=Indel_haplotype))+
  geom_jitter(aes(y=TPM, color=Indel_haplotype),width = 0.2, alpha=0.5)+
  geom_boxplot(aes(y=TPM), alpha=0.75)+
  #geom_text(aes(y=TPM,label=Genome),size=3,alpha=0.6,position = "jitter")+
  facet_wrap(~Tissue, scales = "free")+
  ggtitle(paste("Zm00001e037501",t$Region_Type[1],sep=" "))

#Checking the data before one-way anova of GL15_promoter (Zm00001e037501)
group_by(t, Indel_haplotype) %>%filter(Tissue=="V11_tip") %>%summarise(count = n(), mean=mean(TPM, na.rm = TRUE), sd=sd(TPM,na.rm = TRUE))
group_by(t, Indel_haplotype,Tissue)  %>%summarise(count = n(), mean=mean(TPM, na.rm = TRUE), sd=sd(TPM,na.rm = TRUE)) %>% print(n=30)
#one-way anova
res.aov<-aov(TPM ~ Indel_haplotype, data = filter(t, Tissue=="V11_tip"))
summary(res.aov)
TukeyHSD(res.aov)




#Investigating why there are so many datapoints
#It's because the TPM has transcript names
t[grep("_T001",t$Geneid),] %>% ggplot(aes(x=Indel_haplotype))+
  geom_jitter(aes(y=TPM, color=Indel_haplotype),width = 0.2, alpha=0.5)+
  geom_boxplot(aes(y=TPM), alpha=0.75)+
  #geom_text(aes(y=TPM,label=Genome),size=3,alpha=0.6,position = "jitter")+
  facet_wrap(~Tissue, scales = "free")+
  ggtitle(paste("Zm00001e037501",t$Region_Type[1],sep=" "))

#Haplotype investigation of phyB1 (Zm00001e001462)
#For some reason MS71 expression data is not in the AllTPM table for this gene
phyB1_promogene<-select(`G.Zm00001e001462 chr1:49879953-49891691`, c(Region,Region_Type,Genome,Flower_Time))
phyB1_promogene<-add_column(phyB1_promogene, Indel_haplotype = c("del3270","del3270","none","del3270-del3272","none",
                                                                 "del3270","del3270","del3270","none","ins3273","none",
                                                                 "del3270","none","del3270","none","ins3273","ins3271","none",
                                                                 "ins3273","none","none","del3270-del3272","ins3273","none",
                                                                 "ins3273"))
t<-inner_join(filter(AllTPM, V3 == "GRMZM2G124532_T01"),phyB1_promogene,by="Genome")
t[,"Indel_haplotype"]<-factor(t$Indel_haplotype, levels = c("none","del3270","del3270-del3272","ins3271","ins3273"))
t$Tissue<-factor(t$Tissue, levels=c("V11_base","V11_middle","V11_tip","V18_ear","V18_tassel","R1_anther"),ordered=TRUE)
ggplot(data = t,aes(x=Indel_haplotype))+
  geom_jitter(aes(y=TPM, color=Indel_haplotype),width = 0.2, alpha=0.5)+
  geom_boxplot(aes(y=TPM), alpha=0.75)+
  #geom_text(aes(y=TPM,label=Genome),size=3,alpha=0.6,position = "jitter")+
  facet_wrap(~Tissue, scales = "free")+
  ggtitle(paste("Zm00001e001462",t$Region_Type[1],sep=" "))

#Haplotype investigation of phyA1 (no_evd_annotation)
phyA1_promogene<-select(`G.no_evd_annotation chr1:275855389-275861912`, c(Region,Region_Type,Genome,Flower_Time))
phyA1_promogene<-add_column(phyA1_promogene, Indel_haplotype = c("del14213","del14213","del14215",
                                                                 "del14213-ins14-del15","none","del14215",
                                                                 "none","ins14214-del15","del14213","none",
                                                                 "del14213","none","none","del14215","none",
                                                                 "del14213","del14213","del14213","none",
                                                                 "del14213-del15","del14213","none","ins14214-del15",
                                                                 "none","none"))
t<-inner_join(filter(AllTPM, V3 == "GRMZM2G157727_T01"),phyA1_promogene,by="Genome")
t[,"Indel_haplotype"]<-factor(t$Indel_haplotype, levels = c("none","del14213","del14215","del14213-del15","ins14214-del15","del14213-ins14-del15"))
t$Tissue<-factor(t$Tissue, levels=c("V11_base","V11_middle","V11_tip","V18_ear","V18_tassel","R1_anther"),ordered=TRUE)
ggplot(data = t,aes(x=Indel_haplotype))+
  geom_jitter(aes(y=TPM, color=Indel_haplotype),width = 0.2, alpha=0.5)+
  geom_boxplot(aes(y=TPM), alpha=0.75)+
  geom_text(aes(y=TPM,label=Genome),size=3,alpha=0.6,position = "jitter")+
  facet_wrap(~Tissue, scales = "free")+
  ggtitle(paste("Zm00001e001462",t$Region_Type[1],sep=" "))

#Haplotype investigation of ZmTOC1 (Zm00001e022997)
#There are no interactions between the two indels of interest, and 1 is only present in 1 line
#There does not appear to be a reason to run this analysis at this time
#RECHECK THIS ONE BECAUSE THERE"S 4 INDELS