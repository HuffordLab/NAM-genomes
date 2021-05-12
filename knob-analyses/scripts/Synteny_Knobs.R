##Synteny 
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)

###ARRAY DATA. This is output from Rep_content_dat_prep.R
setwd("~/Desktop/Dawe/Summary")
rep_coords<- read.table("NAM_array_coords.tsv", header = T) ##All Array data from the Rep_content_dat_prep.R script
rep_coords<- rep_coords[rep_coords$V1 != "B73_AB10",] ##Removed Ab10

###SYNTOLOG DATA -- Sorghum genes
setwd("~/Desktop/Dawe/Sorghum_final_annotations_bedfiles") ## where all the Sorghum syntolog files are for me 
path="~/Desktop/Dawe/Sorghum_final_annotations_bedfiles"
beds<- dir(path, pattern=".bed", recursive=TRUE, full.names=TRUE) ##all the summary files
beds_line<- str_split( beds, pattern="-", simplify=T)[,2]
beds_line[12]<-  "IL14H" # rename to match the NAM_array_coords.tsv file
beds_line[19]<-  "MS71"
beds_line[23]<-  "Oh7b" 

###CYTOLOGICAL DATA
setwd("~/Desktop/Dawe/Summary") 
cyt<- read.csv("NAM_array_coords_annotation_cyt_search.csv")

#GENOME CHROMOSOME LENS
#these are just chromosome lengths from the psueodomolecules

path="~/Desktop/Dawe/genomes"
lens<- dir(path, pattern="fasta.len", recursive=TRUE, full.names=TRUE) ##all the summary files
lens<- lens[2:length(lens)] ##remove Ab10 dat
lens<- c( lens[1], lens [3:length(lens)])
B73_lens<- read.table(lens[1])
lens_lines<- as.data.frame(matrix( nrow=length(lens), ncol=11))
lens_nam<- str_split( lens, pattern="/", simplify=T)[,7]
lens_lines$V1<-str_split( lens_nam, pattern="[.]", simplify=T)[,1]

for( i in 1:length(lens)){
  fil<- read.table( lens[i])
  lens_lines[i,2:11]<- fil[1:10,2]
}

##Len files here are just the lengths of the psuedomolecules for plotting relative position
#in a directory containing all the assemblies as fastas, do the following
#ls *.fasta > pseudos

#while read f
#do
#cat $f | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $f.len
#done < pseudos


#just some useful chr
chr_num<- chr_num<- paste( "chr", 1:10, sep="")
rep_coords_chr<- rep_coords[rep_coords$Chr %in% chr_num,]

## Using their respective Sorghum gene bed files, I will get the up and down
# sorghum genes for each array for each line


############################################
### Main figure version-- synteny to significant B73 knobs!
rep_coords_all_synmap<- as.data.frame(matrix(ncol=ncol(rep_coords_chr)+4 , nrow=0))

for(j in 1:length(beds_line)){
  sub_elems<-rep_coords_chr[rep_coords_chr$V1 == beds_line[j] & rep_coords_chr$Chr %in% chr_num ,]
  sub_elems$up<-NA #gene from own line
  sub_elems$down<-NA
  sub_elems$up_Sor<-NA #Sorghum ortholog (connects all the lines together)
  sub_elems$down_Sor<-NA
  
  bed<- read.table(beds[j])
  bed$V7<- paste( bed$V5, bed$V6, sep="_")
  gff<- bed
  
  for(i in 1:nrow(sub_elems)){
    tar<- sub_elems[i,]
    gff_low<-  gff[  gff$V2 < tar$Start & as.character(gff$V1) == as.character(tar$Chr),]
    gff_low_or<- gff_low[order(-gff_low$V2), ]
    sub_elems$up[i]<-as.character(gff_low_or$V4[1])
    sub_elems$up_Sor[i]<-as.character(gff_low_or$V7[1])
    
    gff_high<-  gff[  gff$V3> tar$End & as.character(gff$V1) == as.character(tar$Chr),]
    gff_high_or<- gff_high[order(gff_high$V3), ]
    sub_elems$down[i]<-as.character(gff_high_or$V4[1])
    sub_elems$down_Sor[i]<-as.character(gff_high_or$V7[1])
  }
  
  rep_coords_all_synmap<- rbind(rep_coords_all_synmap,  sub_elems)
  
}

#only the knobs
rep_coords_all_synmap_knob<-rep_coords_all_synmap[  rep_coords_all_synmap$type %in% c("TR-1", "knob180"),]

#########
##For our purposes, we're going to only want knobs that are syntenic to a knob present
# in B73 that is at least 100kb. So let's do that here

rep_coords_all_synmap_knob$uni_ID<- paste(rep_coords_all_synmap_knob$Chr, rep_coords_all_synmap_knob$up_Sor,  rep_coords_all_synmap_knob$down_Sor, sep="_")
B73_filt<- rep_coords_all_synmap_knob[rep_coords_all_synmap_knob$V1 == "B73" & rep_coords_all_synmap_knob$Array.size..bp. >= 100000, ]
B73_ids<- unique(B73_filt$uni_ID)
##note- there are 3 arrays in in B73 that share sytenic spots  (22 rows of B73_filt, 19 B73_ids)

rep_coords_all_synmap_knob_B73_filt<- rep_coords_all_synmap_knob[rep_coords_all_synmap_knob$uni_ID %in% B73_ids, ]
not_B73<- rep_coords_all_synmap_knob_B73_filt[rep_coords_all_synmap_knob_B73_filt$V1 != "B73",] 
is_B73<- rep_coords_all_synmap_knob_B73_filt[rep_coords_all_synmap_knob_B73_filt$V1 == "B73",]
##here, I also make a B73 subset. This is because there are 2 extra arrays in this data versus
# the B73_filt. That tells me there are 2 arrays that share the Sorghum gene synteny with larger 
# arrays but are less than 100kb. (24 rows for is_B73, 22 for B73_filt)

not_sub_pos<- as.data.frame(matrix(nrow=0, ncol=ncol(not_B73)+2))
is_B73$comp_dist<-NA
for( i in 1:length(B73_ids)){
  not_sub<- not_B73[not_B73$uni_ID==B73_ids[i],]
  is_sub<- is_B73[is_B73$uni_ID==B73_ids[i],]
  
  if( nrow(not_sub) > 0  ){
    

  for( k in 1:length(unique(not_sub$V1))){
    not_sub_lin<- not_sub[not_sub$V1 == unique(not_sub$V1)[k],]

  if( nrow(not_sub_lin) == nrow(is_sub)){
    ##if B73 and the line of interest have the same number of arrays at the syntenic spot of interest,
    # then we assume they are sytenic in the same order
    
      not_sub_lin$pos_start<- NA
      not_sub_lin$pos_start<- is_sub$Start
      not_sub_lin$pos_end<-NA
      not_sub_lin$pos_end<- is_sub$End
      
      not_sub_pos<- rbind(not_sub_pos, not_sub_lin)
  } else {
    
    ##otherwise, this section determines the similar array by the absolutel value of the 
    #differences in the start coordinates
    
    for( j in 1:nrow(not_sub_lin)){
      tar<- not_sub_lin[j,]
    
      is_sub$comp_dist<- abs(is_sub$Start - tar$Start[1])
      is_sub_or<- is_sub[order(is_sub$comp_dist),] 
    
      tar$pos_start<- is_sub_or$Start[1]
      tar$pos_end<- is_sub_or$End[1]
    
      not_sub_pos<- rbind(not_sub_pos, tar)

      }
  }
  
  }
  }
  
}

sub_dat_not<- select(not_sub_pos, c("V1", "Chr", "pos_start", "pos_end" , "type", "Array.size..bp." ))
sub_dat_is<- select(is_B73, c("V1", "Chr", "Start", "End" , "type", "Array.size..bp." ))

nams<- c("Line", "Chr","Start", "End", "type", "original_size")
colnames(sub_dat_not)<- nams
colnames(sub_dat_is)<- nams

new_dat<- rbind(sub_dat_is,sub_dat_not)


new_dat_or<- new_dat[order(new_dat$Chr, new_dat$Start),]

write.csv(new_dat_or, "knob_synteny_chr_B73_100kb_v2.csv", row.names = F)

#Finished data for main Figure!
############################################

### Supplement figure version-- pansynteny for all lines!
#To compare to cytological elements

## cut the files down to a common set of genes that are syntenic in all
B73<- read.table(beds[1])
B73$V7<- paste( B73$V5, B73$V6, sep="_")
out<- B73$V7

for( i in 2:length( beds)){
  tab<- read.table(beds[i])
  tab$V7 <- paste( tab$V5, tab$V6, sep="_" )
  out<- intersect( out, tab$V7)
}

## length(out) = 18539
##out is a common set of genes for all 
rep_coords_all_synmap<- as.data.frame(matrix(ncol=ncol(rep_coords_chr)+4 , nrow=0))

for(j in 1:length(beds_line)){
  sub_elems<-rep_coords_chr[rep_coords_chr$V1 == beds_line[j] & rep_coords_chr$Chr %in% chr_num ,]
  sub_elems$up<-NA #gene from own line
  sub_elems$down<-NA
  sub_elems$up_coord<-NA #Sorghum ortholog (connects all the lines together)
  sub_elems$down_coord<-NA
  
  bed<- read.table(beds[j])
  bed$V7<- paste( bed$V5, bed$V6, sep="_")
  gff<- bed[bed$V7 %in% out, ]
  
  for(i in 1:nrow(sub_elems)){
    tar<- sub_elems[i,]
    gff_low<-  gff[  gff$V2 < tar$Start & as.character(gff$V1) == as.character(tar$Chr),]
    gff_low_or<- gff_low[order(-gff_low$V2), ]
    sub_elems$up[i]<-as.character(gff_low_or$V7[1])
    sub_elems$up_coord[i]<-as.character(gff_low_or$V3[1]) ##end coord from up gene
    
    gff_high<-  gff[  gff$V3> tar$End & as.character(gff$V1) == as.character(tar$Chr),]
    gff_high_or<- gff_high[order(gff_high$V3), ]
    sub_elems$down[i]<-as.character(gff_high_or$V7[1])
    sub_elems$down_coord[i]<-as.character(gff_high_or$V2[1]) ##start coord from down gene
  }
  
  rep_coords_all_synmap<- rbind(rep_coords_all_synmap,  sub_elems)
  
}

#add relative coordinates from own line
rep_coords_all_synmap$up_coord_rel<-NA
rep_coords_all_synmap$down_coord_rel<-NA

#add the relative positions from B73
rep_coords_all_synmap$up_coord_B73<- NA
rep_coords_all_synmap$down_coord_B73<- NA
rep_coords_all_synmap$up_coord_relB73<- NA
rep_coords_all_synmap$down_coord_relB73<- NA

for( i in 1:nrow( rep_coords_all_synmap)) {
  
  tar<- lens_lines[lens_lines$V1 == rep_coords_all_synmap$V1[i], ]
  chr<- as.numeric( str_split(rep_coords_all_synmap$Chr[i], pattern="r", simplify=T)[,2]) 
  
  #fix empty coords
  if( is.na( rep_coords_all_synmap$up_coord[i])){  rep_coords_all_synmap$up_coord[i] <- 1} # if there is no start coord, enter 1
  if( is.na( rep_coords_all_synmap$down_coord[i])){ # if there is no end coord, enter chromosome len
    rep_coords_all_synmap$down_coord[i] <- tar[1, chr+1]}
  
  #relative position ( up_coord_rel = up_coord/chromosome len)
  if(rep_coords_all_synmap$up_coord[i]>1){ rep_coords_all_synmap$up_coord_rel[i] <- as.numeric(rep_coords_all_synmap$up_coord[i] )/tar[1,chr+1] 
  }else{ rep_coords_all_synmap$up_coord_rel[i] <-0} #if the up coord is >1 (not the very start of chromosome), then set as zero
  rep_coords_all_synmap$down_coord_rel[i]  <- as.numeric(rep_coords_all_synmap$down_coord[i]) /tar[1,chr+1]

  if( !is.na(rep_coords_all_synmap$up[i]) & rep_coords_all_synmap$up_coord[i] > 1){
    B73_gene_up<- B73[B73$V7 == rep_coords_all_synmap$up[i], ]
    rep_coords_all_synmap$up_coord_B73[i]<- B73_gene_up$V3
    rep_coords_all_synmap$up_coord_relB73[i]<- B73_gene_up$V3/lens_lines[1,chr+1]
  } else {
    rep_coords_all_synmap$up_coord_B73[i]<-1
    rep_coords_all_synmap$up_coord_relB73[i]<-0} 
  #relative coord from B73. Pull corresponding information from B73 Sorghum ortho genes
  
  if( !is.na(rep_coords_all_synmap$down[i])){
    B73_gene_down<- B73[B73$V7 == rep_coords_all_synmap$down[i], ]
    rep_coords_all_synmap$down_coord_B73[i]<- B73_gene_down$V2
    rep_coords_all_synmap$down_coord_relB73[i]<- B73_gene_down$V3/lens_lines[1,chr+1]
  }else { 
    rep_coords_all_synmap$down_coord_relB73[i]<-1
    rep_coords_all_synmap$down_coord_B73[i]<-lens_lines[1,chr+1]}
  #same for down
}


#write.table(rep_coords_all_synmap, "Synt_v2_all_arrays.tsv", sep="\t", quote=F)
#all the synteny basic information

#get the basic column names in agreement for easier appending-- here, the line column seems to be the only not in agreement
colnames( rep_coords_all_synmap)<- c("NAM_line", colnames(rep_coords_all_synmap)[2:length(colnames(rep_coords_all_synmap))])
###add in that cytological data

all_dat<-merge(cyt,rep_coords_all_synmap,by=c("NAM_line", "Chr", "Start", "End","type","structure"), all=TRUE) 
##the cyt data has 453 
##the rep_coords_all_synmap data has 2499 becuase it includes all array types
##all_dat will reduce down to 1660, which will be knobs on chromosomes.
all_dat[is.na(all_dat$down),"down"]<- paste( "End", all_dat[is.na(all_dat$down),2], sep="_")
all_dat[is.na(all_dat$up),"up"]<- paste( "Start", all_dat[is.na(all_dat$up),2], sep="_")

#just subset down the chromosomes for the rest -- 1660 rows
all_dat<- all_dat[all_dat$type %in% c("knob180", "TR-1"),] 

######
# okay, now let's fix the issue spots
# with this high of density, there are some SV apparent
trouble<- list()
trouble_index<-1

all_dat_edit_v1<- as.data.frame( matrix( nrow=0, ncol=ncol(all_dat)+1))
# so I can recombine without missing something later, I'll make a new df with the good guys

all_dat$group<- paste(all_dat$up, all_dat$down, sep="_")

for( i in 1:length( unique(all_dat$up))){
  set<- all_dat[all_dat$up== unique(all_dat$up)[i],] #groups are defined by their start coordinate 
  
  check<- set$group == set$group[1] # check if all in group share that ID (same up/down gene coords)
  
  if( sum(check) != length(check )){  #the sum function counts the number of TRUE's when applied to a logical vector 
    # if not, add that group to the trouble list!
    trouble[[trouble_index]] <- set
    trouble_index<- trouble_index +1
    } else { all_dat_edit_v1<- rbind( all_dat_edit_v1, set)}

}

all_dat_edit<- as.data.frame( matrix( nrow=0, ncol=ncol(all_dat_edit_v1)))
# reset the good guys for the other set of issues

for( i in 1:length( unique(all_dat$down))){
  set<- all_dat_edit_v1[all_dat_edit_v1$down== unique(all_dat$down)[i],] #groups are defined by their start coordinate 
  
  check<- set$group == set$group[1] # check if all in group share that ID (same up/down gene coords)
  
  if( sum(check) != length(check )){  #the sum function counts the number of TRUE's when applied to a logical vector 
    # if not, add that group to the trouble list!
    trouble[[trouble_index]] <- set
    trouble_index<- trouble_index +1
  } else { all_dat_edit<- rbind( all_dat_edit, set)}
  
}

#make sure that there isn't an issue-- there are the expected number of lines
lin_check<- 0
for( i in 1: length(trouble)){
  lin_check<- lin_check + nrow( trouble[[i]])
}

lin_check #384 and the all_dat_edit file is 1276 = 1660 rows, as expected


all_dat_edit$up_fix<-all_dat_edit$up
all_dat_edit$down_fix<-all_dat_edit$down
all_dat_edit$up_coord_B73_fix<-all_dat_edit$up_coord_B73
all_dat_edit$down_coord_B73_fix<-all_dat_edit$down_coord_B73
all_dat_edit$up_coord_relB73_fix<-all_dat_edit$up_coord_relB73
all_dat_edit$down_coord_relB73_fix<-all_dat_edit$down_coord_relB73

## modify the issue coordinates. In this case, most of the trouble spots have one B73 element, 
# which is the first row
for( i in 1:length( trouble)){
  dat<- trouble[[i]]
  if(dat$NAM_line[1]  =="B73" | str_split(dat$up[1], pattern="_", simplify=T )[,1] == "Start"){
  dat$up_fix<- dat$up[1]
  dat$down_fix<- dat$down[1]
  dat$up_coord_B73_fix<- dat$up_coord_B73[1]
  dat$down_coord_B73_fix<- dat$down_coord_B73[1]
  dat$up_coord_relB73_fix<- dat$up_coord_relB73[1]
  dat$down_coord_relB73_fix<- dat$down_coord_relB73[1]
  }else{
    tar_up<- B73[B73$V7 %in% dat$up[1], ]
    new<- B73[B73$V1 == tar_up$V1 & B73$V2 >= tar_up$V2, ]
    new_or<- new[ order( new$V2),]
    tar_down<-new_or[2,]
      
    dat$up_fix<- new_or$V7[1]
    dat$down_fix<- tar_down$V7[1]
    dat$up_coord_B73_fix<- new_or$V3[1]
    dat$down_coord_B73_fix<-  tar_down$V2[1]
    dat$up_coord_relB73_fix<- new_or$V3[1]/B73_lens[B73_lens$V1 == new_or$V1[1],2]
    dat$down_coord_relB73_fix<-  tar_down$V2[1]/B73_lens[B73_lens$V1 == new_or$V1[1],2]
  }
  all_dat_edit<- rbind( all_dat_edit, dat)
  
}

##plotting just the cytological groups
cyt_dat<- all_dat_edit[all_dat_edit$Cyto_vis == "YES", "group"]
cyt_groups<- unique(cyt_dat)
new_dat_or_cyt<- all_dat_edit[all_dat_edit$group %in% cyt_groups,]

new_dat_or_cyt$colour<-NA
for( i in 1:nrow(new_dat_or_cyt)){
  if( new_dat_or_cyt$Cyto_vis[i] %in% "YES"){
    new_dat_or_cyt$colour[i] <- new_dat_or_cyt$type[i]
  } else{ new_dat_or_cyt$colour[i] <- "Not Visible"}
}


new_dat_or_cyt$Chr<- as.factor(new_dat_or_cyt$Chr)
new_dat_or_cyt$Chr<- factor(new_dat_or_cyt$Chr, levels( new_dat_or_cyt$Chr)[c( 1,3:10,2)])

new_dat_or_cyt$NAM_line<- as.factor(new_dat_or_cyt$NAM_line)
new_dat_or_cyt$NAM_line<- factor(new_dat_or_cyt$NAM_line, levels = c("B73", "B97",  "Ky21", "M162W", "MS71", "Oh43", "Oh7b",
                           "M37W", "Mo18W" , "Tx303", "HP301","P39", "IL14H", "CML52", "CML69",
                           "CML103" , "CML228", "CML247", "CML277", "CML322", "CML333", "Ki3",
                           "Ki11" , "NC350" , "NC358", "Tzi8"))
new_dat_or_cyt$NAM_line<- factor(new_dat_or_cyt$NAM_line, rev(levels(new_dat_or_cyt$NAM_line)))

x_cols <- c(rep(c("goldenrod1"),1),rep(c("gray47"),3),rep(c("royalblue"),6),rep(c("orchid"),1),rep(c("orangered"),2),rep(c("limegreen"),13))
x_cols2<- rev(x_cols)

p7 <- ggplot(new_dat_or_cyt, aes( x=up_coord_relB73_fix, y=NAM_line)) +
  geom_line(data=new_dat_or_cyt, aes(x=up_coord_relB73_fix, y=NAM_line, group=down_coord_relB73_fix))  +theme_classic() + facet_wrap(~Chr, nrow=2, strip.position="top") +
  
  geom_point(data=new_dat_or_cyt, aes(x=up_coord_relB73_fix, y=NAM_line, group=down_coord_relB73_fix ,colour=colour, size= Array.size..bp.)) + 
  scale_colour_manual(values=c('cadetblue2',  'black', 'deeppink3')) +
  labs( x="Relative Position in B73",y="Line")  +
  theme(axis.text.y=element_text(size=15,hjust = 1, colour= x_cols2),text = element_text(size=15), panel.background = element_blank(), legend.position = "none")
p7

ggsave( plot=p7 , "Supplement_cyto_knobs_2.png", dpi=600, width=15, height=10)
