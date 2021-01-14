###Repeat Array data prep


#########This section of the code is to reorganize the Repeat_content_sum.csv files 

#All Rep content files downloaded from the project folder to my local folder
setwd("~/Desktop/Dawe/Summary") ##where all the Repeat_content_sum.csv files are held
path = "~/Desktop/Dawe/Summary"

library("dplyr")
library("stringr")
library("stats")

##data below are the output from the array identification python script. 
#Since they all have the same name, Repeat_content_sum.csv, they were
#copied to my local computer for running this script with an added 
#suffix of "_<line name>"

#DATA FROM ARRAYS
data_files<- dir(path, pattern="Repeat_content_sum.csv_", recursive=TRUE, full.names=TRUE) ##all the summary files

#linux code 
# #copy all the needed Repeat_content_sum.csv files to local directory

#cd $DIR #directory for the data
#mkdir rename_Rep
#ls -ltrh */Repeat_conte* > fils
#for f in */Repeat_conte* 
#  do
#cp ./$f ./rename_Rep/"${f#*/}"_"${f%%/*}"
#done

##copy all those Repeat_content_sum.csv_ files to local folder. In my case, ~/Desktop/Dawe/Summary



elem_names<- c( "AF013103.1", "CentC", "CL569186.1", "knob180", "TR-1")


All <- lapply(data_files,function(i){
  read.csv(i, header=TRUE)
}) #Read all the Repeat_content_sum files in here, put in a list object called 'All'


n2<-str_split(data_files, pattern="/", simplify=T)
NAM_names<- n2[,7]
#make a vector of the line names in the same order as they are in the All list. 


chrs<- 1:10
chr_num<- paste( "chr", chrs, sep="") #this will come in handy later when you want to look just at chromosome elements


#convert to summary data table with mixing information that combines and simplifies information from the original data structure

rep_mixing<- function(k) {
  
  Rep_content<-k
  
  add_col<-as.data.frame(matrix(ncol=5,nrow=nrow(Rep_content))) 
  
  colnames(add_col)<- c("AF013103.1", "CentC", "CL569186.1",   
                        "knob180", "TR-1"  )
  
  for( i in 1:nrow(Rep_content)){
    dat<-Rep_content[i, 12:14]
    colnames(dat)<- c("element", "orientation", "bp")
    add_elem<- (ncol(Rep_content)-14)/3
    
    for(y in 1:add_elem){
      n<- (14 + y*3)
      if(is.na(Rep_content[i, n]) == FALSE){ 
        new<- Rep_content[i, (n-2):n]
        colnames(new)<- c("element", "orientation", "bp")
        dat<- rbind( dat, new)}
    }
    
    
    for( z in 1:5){
      name<- elem_names[z]
      elem<- dat[dat$element == name ,]
      length<- sum(elem$bp)
      length -> add_col[i,z]
      
    }
  }
  
  total_rep_length<-Rep_content$Repeat.content..bp.
  percent_dat<- round(add_col/total_rep_length,3)
  
  Rep_select<- select(Rep_content, c("Chr", "Group", "Start", "End", "Array.size..bp.", 
                                     "Repeat.content..bp.", "Repeat.Percentage", "X100N", "Ngap.Percentage"))
  
  
  total_select<- cbind(Rep_select, percent_dat)
  
  total_sel_filt<- total_select[total_select$Repeat.Percentage >= 0.10, ]
  
  return(total_sel_filt)
}

el<-lapply(All,rep_mixing)


#### if you want to group all the el into on big df (all elements together, can be simplified and filtered for just the elements you want)
#######
el_named<- list()

for( i in 1:length(el)){
  nam<- as.data.frame(matrix(nrow=nrow(el[[i]]), ncol=1))
  nam[,1]<- NAM_names[i]
  colnames(nam)<- "NAM_line"
  el_named[[i]]<- cbind(nam, el[[i]])
}

Rep_arrays_data<- el_named[[1]]
for( i in 2:length(el_named)){
  Rep_arrays_data<- rbind( Rep_arrays_data, el_named[[i]])
}

##To separate combined file back to a list object, if needed
NAM_names_unique<- unique( Rep_arrays_data$NAM_line)

el_named_relist<- list()
for( i in length(NAM_names_unique)){
  dat<- Rep_arrays_data[Rep_arrays_data$NAM_line == nam,]
  dat<- select( dat, -c("NAM_line"))
  el_named_relist[[i]]<- dat 
}
########

#I save all of el as one big df for ease. 
write.table(Rep_arrays_data, "Rep_arrays_data.tsv", sep="\t",  col.names = T, row.names = F, quote = F)
#all arrays, but they aren't labeled y major repeat yet


#######
#filtering is for all repeat arrays that are 100kb with at least 10% repeat content


rep_filt<- function(j) {
  Prep_dat<-j
  sub<- Prep_dat[Prep_dat$Array.size..bp. >= 1, ] #100kb+
  sub2<- sub[sub$Repeat.Percentage>=.10 & sub$Array.size..bp.> 10000, ] #knob180 and TR-1 arrays
  sub2$type<-NA
  sub2$structure<-NA
  for( i in 1:nrow(sub2)){
    if( (sub2$knob180[i]  + sub2$'TR-1'[i]) > .50 ){
        if( sub2$knob180[i] > sub2$'TR-1'[i]) { sub2$type[i]<- "knob180" }
          else {sub2$type[i]<- "TR-1" } }
    else if( sub2$CentC[i] > .50 ){ sub2$type[i]<- "CentC" }
    else if( sub2$AF013103.1[i] > .50 ){ sub2$type[i]<- "AF013103.1" }
    else if( sub2$CL569186.1[i] > .50 ){ sub2$type[i]<- "CL569186.1" }
    else sub2$type[i]<- "NA" 
    
    values<-sub2[i,10:14]
    ifelse( sum(values > 0) >=2, sub2$structure[i]<- "mixed", sub2$structure[i]<- "single")
  }
  
  
  
  return(sub2)
}


filt_el<- lapply(el,rep_filt) #filtered and labeled


#make sure all have good labels

sum_NA<- as.data.frame(matrix(nrow=length(filt_el), ncol=1))

for( i in 1:length(filt_el)){
  dat<- filt_el[[i]]
  s<- sum(dat$type == "NA")
  sum_NA[i,1]<-s
}

##combine into one again

lab_filt_el<- list()


for( i in 1:length(NAM_names)){
  dat<- filt_el[[i]]
  name<- NAM_names[i]
  
  line<- as.data.frame(matrix(nrow=nrow(dat), ncol=1))
  line[,1]<- name
  lab_dat<- cbind(line, dat)
  
  lab_filt_el[[i]]<-lab_dat[lab_dat$Array.size..bp.>= 10000 , ]
}

comb_coords<- lab_filt_el[[1]] #this will be the summary and coordinate information for all knobs annotated 


for(i in 2:length(NAM_names)){
  comb_coords<- rbind(comb_coords, lab_filt_el[[i]])
}

comb_coords$V1<- str_split(comb_coords$V1, pattern=".csv_", simplify = T)[,2]

write.table(comb_coords, "NAM_array_coords.tsv", quote = F, row.names = F, col.names = T) 
#This is the base file for downstream analyses

## note that this data also includes CentC arrays as identified, but are not necessarily correlated to the function centromere array data 
## prepared by Jianing Liu