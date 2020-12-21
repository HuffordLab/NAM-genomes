library(dplyr)
library(plyr)
library(magrittr)
# set dir for the whole work,where the files are located 
setwd("~/Desktop/NAM_26_genomes/nucmer_1000_10/base_files/")
#reading canonical sequences from all name genomes. These will serve as the base pan gene for each lines

B73c <- read.csv("B73_canonical_transcript.txt", sep = "\t")
B97c <- read.csv("B97_canonical_transcript.txt", sep = "\t")
Ky21c <- read.csv("Ky21_canonical_transcript.txt", sep = "\t")
M162Wc <- read.csv("M162W_canonical_transcript.txt", sep = "\t")
Ms71c <- read.csv("Ms71_canonical_transcript.txt", sep = "\t")
Oh7Bc <- read.csv("Oh7B_canonical_transcript.txt", sep = "\t")
Oh43c <- read.csv("Oh43_canonical_transcript.txt", sep = "\t")
M37Wc <- read.csv("M37W_canonical_transcript.txt", sep = "\t")
Mo18Wc <- read.csv("Mo18W_canonical_transcript.txt", sep = "\t")
Tx303c <- read.csv("Tx303_canonical_transcript.txt", sep = "\t")
HP301c <- read.csv("HP301_canonical_transcript.txt", sep = "\t")
Il14Hc <- read.csv("Il14H_canonical_transcript.txt", sep = "\t")
P39c <- read.csv("P39_canonical_transcript.txt", sep = "\t")
CML52c <- read.csv("CML52_canonical_transcript.txt", sep = "\t")
CML69c <- read.csv("CML69_canonical_transcript.txt", sep = "\t")
CML103c <- read.csv("CML103_canonical_transcript.txt", sep = "\t")
CML228c <- read.csv("CML228_canonical_transcript.txt", sep = "\t")
CML247c <- read.csv("CML247_canonical_transcript.txt", sep = "\t")
CML277c <- read.csv("CML277_canonical_transcript.txt", sep = "\t")
CML322c <- read.csv("CML322_canonical_transcript.txt", sep = "\t")
CML333c <- read.csv("CML333_canonical_transcript.txt", sep = "\t")
Ki3c <- read.csv("Ki3_canonical_transcript.txt", sep = "\t")
Ki11c <- read.csv("Ki11_canonical_transcript.txt", sep = "\t")
NC350c <- read.csv("NC350_canonical_transcript.txt", sep = "\t")
NC358c <- read.csv("NC358_canonical_transcript.txt", sep = "\t")
Tzi8c <- read.csv("Tzi8_canonical_transcript.txt", sep = "\t")

# prepare two column, self duplicated gene list as left joint input for individual NAM line 
B73_self <- read.csv("B73.txt", sep = "\t")
P39_self <- read.csv("P39.txt", sep = "\t")
Ki11_self <- read.csv("Ki11.txt", sep = "\t")
NC350_self <- read.csv("NC350.txt", sep = "\t")
Oh43_self <- read.csv("Oh43.txt", sep = "\t")
Oh7B_self <- read.csv("Oh7B.txt", sep = "\t")
Ky21_self <- read.csv("Ky21.txt", sep = "\t")
Ms71_self <- read.csv("Ms71.txt", sep = "\t")
M162W_self <- read.csv("M162W.txt", sep = "\t")
Il14H_self <- read.csv("Il14H.txt", sep = "\t")
HP301_self <- read.csv("HP301.txt", sep = "\t")
CML52_self <- read.csv("CML52.txt", sep = "\t")
CML69_self <- read.csv("CML69.txt", sep = "\t")
CML228_self <- read.csv("CML228.txt", sep = "\t")
CML247_self <- read.csv("CML247.txt", sep = "\t")
CML277_self <- read.csv("CML277.txt", sep = "\t")
CML322_self <- read.csv("CML322.txt", sep = "\t")
Mo18W_self <- read.csv("Mo18W.txt", sep = "\t")
Tzi8_self <- read.csv("Tzi8.txt", sep = "\t")
M37W_self <- read.csv("M37W.txt", sep = "\t")
NC358_self <- read.csv("NC358.txt", sep = "\t")
B97_self <- read.csv("B97.txt", sep = "\t")
Tx303_self <- read.csv("Tx303.txt", sep = "\t")
Ki3_self <- read.csv("Ki3.txt", sep = "\t")
CML103_self <- read.csv("CML103.txt", sep = "\t")
CML333_self <- read.csv("CML333.txt", sep = "\t")


# adding the first NAM genome: B73 
## pan 1 
Pan_reference <- read.csv("B73_canonical_transcript.txt")
B73 <- B73_self
B97 <- read.csv("B73_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("B73_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("B73_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("B73_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("B73_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("B73_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("B73_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("B73_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("B73_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("B73_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("B73_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("B73_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("B73_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("B73_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("B73_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("B73_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("B73_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("B73_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("B73_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("B73_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("B73_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("B73_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("B73_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("B73_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("B73_Tzi8_all_pair.txt", sep = "\t")
Filtered_B73_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_B73_pan,file="Filtered_B73_pan.csv")

# prepare files for the second genome Tzi8 
# get a list of Tzi8 gene present in Filtered_B73_pan
# get Tzi8 present list 
Tzi8_observers <- as.matrix(unique(unlist(strsplit(as.character(Filtered_B73_pan$Tzi8), ";"))))
Tzi8_pan_id <- as.matrix(Tzi8c[ ! Tzi8c$Query_gene %in% Tzi8_observers[,1], ] )
colnames(Tzi8_pan_id) <- "Query_gene"
write.csv(Tzi8_pan_id, file = "Pan_Tzi8_id.txt")

# adding the 2nd NAM genome: Tzi8 
## pan2 
Pan_reference = read.csv("Pan_Tzi8_id.txt")
Pan_reference[,2]
Tzi8 <- Tzi8_self
B73 <- read.csv("Tzi8_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Tzi8_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Tzi8_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Tzi8_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Tzi8_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Tzi8_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Tzi8_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Tzi8_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Tzi8_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Tzi8_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Tzi8_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Tzi8_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Tzi8_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Tzi8_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Tzi8_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Tzi8_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Tzi8_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Tzi8_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Tzi8_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Tzi8_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Tzi8_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Tzi8_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Tzi8_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Tzi8_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Tzi8_NC358_all_pair.txt", sep = "\t")
Filtered_Tzi8_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Tzi8_pan[-1],file="Filtered_Tzi8_pan.csv")
# prepare files for the 3rd genome Ky21 
# joint dataframe from previous two steps, watch out #variable change 
pan2 <- rbind(Filtered_B73_pan,Filtered_Tzi8_pan[-1])
post_pan2_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan2$Ky21), ";"))))
Ky21_pan_id <- as.matrix(Ky21c[ ! Ky21c$Query_gene %in% post_pan2_remove_list[,1], ] )
colnames(Ky21_pan_id) <- "Query_gene"
write.csv(Ky21_pan_id, file = "Pan_Ky21_id.txt")

# adding the 3rd NAM genome: Ky21 
# pan3 
Pan_reference <- read.csv("Pan_Ky21_id.txt")
Pan_reference[,2]
Ky21 <- Ky21_self
B73 <- read.csv("Ky21_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Ky21_B97_all_pair.txt", sep = "\t")
M162W <- read.csv("Ky21_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Ky21_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Ky21_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Ky21_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Ky21_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Ky21_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Ky21_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Ky21_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Ky21_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Ky21_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Ky21_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Ky21_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Ky21_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Ky21_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Ky21_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Ky21_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Ky21_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Ky21_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Ky21_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Ky21_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Ky21_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Ky21_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Ky21_Tzi8_all_pair.txt", sep = "\t")
Filtered_Ky21_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Ky21_pan, file="Filtered_Ky21_pan.csv")
# prepare files for the 4th genome M162W 
# joint dataframe from previous two steps, watch out #variable change 
pan3 <- rbind(pan2,Filtered_Ky21_pan[,-1])
post_pan3_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan3$M162W), ";"))))
M162W_pan_id <- as.matrix(M162Wc[ ! M162Wc$Query_gene %in% post_pan3_remove_list[,1], ] )
colnames(M162W_pan_id) <- "Query_gene"
write.csv(M162W_pan_id, file = "Pan_M162W_id.txt")

# adding the 4th NAM genome: M162W 
# pan4
Pan_reference <- read.csv("Pan_M162W_id.txt")
Pan_reference[,2]
M162W <- M162W_self
B73 <- read.csv("M162W_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("M162W_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("M162W_Ky21_all_pair.txt", sep = "\t")
Ms71 <- read.csv("M162W_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("M162W_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("M162W_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("M162W_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("M162W_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("M162W_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("M162W_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("M162W_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("M162W_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("M162W_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("M162W_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("M162W_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("M162W_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("M162W_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("M162W_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("M162W_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("M162W_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("M162W_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("M162W_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("M162W_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("M162W_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("M162W_Tzi8_all_pair.txt", sep = "\t")
Filtered_M162W_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_M162W_pan, file="Filtered_M162W_pan.csv")
# prepare files for the 5th genome Ms71 
# joint dataframe from previous two steps, watch out #variable change 
pan4 <- rbind(pan3,Filtered_M162W_pan[,-1])
post_pan4_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan4$Ms71), ";"))))
Ms71_pan_id <- as.matrix(Ms71c[ ! Ms71c$Query_gene %in% post_pan4_remove_list[,1], ] )
colnames(Ms71_pan_id) <- "Query_gene"
write.csv(Ms71_pan_id, file = "Pan_Ms71_id.txt")


# adding the 5th NAM genome: Ms71 
# pan5
Pan_reference <- read.csv("Pan_Ms71_id.txt")
Pan_reference[,2]
Ms71 <- Ms71_self
B73 <- read.csv("Ms71_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Ms71_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Ms71_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Ms71_M162W_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Ms71_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Ms71_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Ms71_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Ms71_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Ms71_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Ms71_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Ms71_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Ms71_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Ms71_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Ms71_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Ms71_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Ms71_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Ms71_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Ms71_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Ms71_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Ms71_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Ms71_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Ms71_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Ms71_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Ms71_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Ms71_Tzi8_all_pair.txt", sep = "\t")
Filtered_Ms71_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Ms71_pan, file="Filtered_Ms71_pan.csv")
# prepare files for the 6th genome Oh7B
# joint dataframe from previous two steps, watch out #variable change 
pan5 <- rbind(pan4,Filtered_Ms71_pan[,-1])
post_pan5_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan5$Oh7B), ";"))))
Oh7B_pan_id <- as.matrix(Oh7Bc[ ! Oh7Bc$Query_gene %in% post_pan5_remove_list[,1],] )
colnames(Oh7B_pan_id) <- "Query_gene"
write.csv(Oh7B_pan_id, file = "Pan_Oh7B_id.txt")


# adding the 6th NAM genome: Oh7B  
# pan 6 
Pan_reference <- read.csv("Pan_Oh7B_id.txt")
Pan_reference[,2]
Oh7B <- Oh7B_self
B73 <- read.csv("Oh7B_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Oh7B_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Oh7B_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Oh7B_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Oh7B_Ms71_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Oh7B_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Oh7B_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Oh7B_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Oh7B_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Oh7B_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Oh7B_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Oh7B_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Oh7B_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Oh7B_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Oh7B_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Oh7B_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Oh7B_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Oh7B_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Oh7B_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Oh7B_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Oh7B_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Oh7B_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Oh7B_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Oh7B_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Oh7B_Tzi8_all_pair.txt", sep = "\t")
Filtered_Oh7B_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Oh7B_pan, file="Filtered_Oh7B_pan.csv")
# prepare files for the 7th genome Oh43
# joint dataframe from previous two steps, watch out #variable change 
pan6 <- rbind(pan5,Filtered_Oh7B_pan[,-1])
post_pan6_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan6$Oh43), ";"))))
Oh43_pan_id <- as.matrix(Oh43c[ ! Oh43c$Query_gene %in% post_pan6_remove_list[,1],] )
colnames(Oh43_pan_id) <- "Query_gene"
write.csv(Oh43_pan_id, file = "Pan_Oh43_id.txt")


# adding the 7th NAM genome: Oh43 
# pan 7 
# adding Oh43 
Pan_reference <- read.csv("Pan_Oh43_id.txt")
Pan_reference[,2]
Oh43 <- Oh43_self
B73 <- read.csv("Oh43_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Oh43_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Oh43_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Oh43_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Oh43_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Oh43_Oh7B_all_pair.txt", sep = "\t")
M37W <- read.csv("Oh43_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Oh43_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Oh43_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Oh43_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Oh43_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Oh43_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Oh43_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Oh43_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Oh43_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Oh43_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Oh43_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Oh43_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Oh43_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Oh43_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Oh43_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Oh43_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Oh43_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Oh43_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Oh43_Tzi8_all_pair.txt", sep = "\t")
Filtered_Oh43_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Oh43_pan, file="Filtered_Oh43_pan.csv")

# prepare files for the 8th genome M37W
# joint dataframe from previous two steps, watch out #variable change 
pan7 <- rbind(pan6,Filtered_Oh43_pan[,-1])
post_pan7_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan7$M37W), ";"))))
M37W_pan_id <- as.matrix(M37Wc[ ! M37Wc$Query_gene %in% post_pan7_remove_list[,1],] )
colnames(M37W_pan_id) <- "Query_gene"
write.csv(M37W_pan_id, file = "M37W_pan_id.txt")


# adding the 8th NAM genome: M37W 
# pan 8 
# add M37W
Pan_reference <- read.csv("M37W_pan_id.txt")
Pan_reference[,2]
M37W <- M37W_self
B73 <- read.csv("M37W_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("M37W_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("M37W_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("M37W_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("M37W_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("M37W_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("M37W_Oh43_all_pair.txt", sep = "\t")
Mo18W <- read.csv("M37W_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("M37W_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("M37W_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("M37W_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("M37W_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("M37W_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("M37W_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("M37W_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("M37W_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("M37W_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("M37W_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("M37W_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("M37W_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("M37W_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("M37W_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("M37W_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("M37W_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("M37W_Tzi8_all_pair.txt", sep = "\t")
Filtered_M37W_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_M37W_pan, file="Filtered_M37W_pan.csv")
# prepare files for the 9th genome Mo18W
# joint dataframe from previous two steps, watch out #variable change 
pan8 <- rbind(pan7,Filtered_M37W_pan[,-1])
post_pan8_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan8$Mo18W), ";"))))
Mo18W_pan_id <- as.matrix(Mo18Wc[ ! Mo18Wc$Query_gene %in% post_pan8_remove_list[,1],] )
colnames(Mo18W_pan_id) <- "Query_gene"
write.csv(Mo18W_pan_id, file = "Mo18W_pan_id.txt")


# adding the 9th NAM genome: Mo18W
# pan 9 
Pan_reference <- read.csv("Mo18W_pan_id.txt")
Pan_reference[,2]

Mo18W <- Mo18W_self
B73 <- read.csv("Mo18W_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Mo18W_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Mo18W_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Mo18W_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Mo18W_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Mo18W_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Mo18W_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Mo18W_M37W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Mo18W_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Mo18W_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Mo18W_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Mo18W_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Mo18W_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Mo18W_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Mo18W_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Mo18W_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Mo18W_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Mo18W_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Mo18W_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Mo18W_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Mo18W_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Mo18W_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Mo18W_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Mo18W_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Mo18W_Tzi8_all_pair.txt", sep = "\t")
Filtered_Mo18W_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Mo18W_pan, file="Filtered_Mo18W_pan.csv")
# prepare files for the 10th genome NC350
# joint dataframe from previous two steps, watch out #variable change 
pan9 <- rbind(pan8,Filtered_Mo18W_pan[,-1])
post_pan9_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan9$NC350), ";"))))
NC350_pan_id <- as.matrix(NC350c[ ! NC350c$Query_gene %in% post_pan9_remove_list[,1],] )
colnames(NC350_pan_id) <- "Query_gene"
write.csv(NC350_pan_id, file = "NC350_pan_id.txt")

# adding the 10th NAM genome: NC350
# pan 10 
# NC350 
Pan_reference <- read.csv("NC350_pan_id.txt")
Pan_reference[,2]

NC350 <- NC350_self
B73 <- read.csv("NC350_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("NC350_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("NC350_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("NC350_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("NC350_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("NC350_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("NC350_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("NC350_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("NC350_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("NC350_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("NC350_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("NC350_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("NC350_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("NC350_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("NC350_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("NC350_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("NC350_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("NC350_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("NC350_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("NC350_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("NC350_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("NC350_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("NC350_Ki11_all_pair.txt", sep = "\t")
NC358 <- read.csv("NC350_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("NC350_Tzi8_all_pair.txt", sep = "\t")
Filtered_NC350_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_NC350_pan, file="Filtered_NC350_pan.csv")
# prepare files for the 11th genome HP301 
# joint dataframe from previous two steps, watch out #variable change 
pan10 <- rbind(pan9,Filtered_NC350_pan[,-1])
post_pan10_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan10$HP301), ";"))))
HP301_pan_id <- as.matrix(HP301c[ ! HP301c$Query_gene %in% post_pan10_remove_list[,1],] )
colnames(HP301_pan_id) <- "Query_gene"
write.csv(HP301_pan_id, file = "HP301_pan_id.txt")

# adding the 11th NAM genome: HP301
# pan 11
Pan_reference <- read.csv("HP301_pan_id.txt")
Pan_reference[,2]
HP301 <- HP301_self
B73 <- read.csv("HP301_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("HP301_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("HP301_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("HP301_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("HP301_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("HP301_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("HP301_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("HP301_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("HP301_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("HP301_Tx303_all_pair.txt", sep = "\t")
Il14H <- read.csv("HP301_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("HP301_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("HP301_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("HP301_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("HP301_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("HP301_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("HP301_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("HP301_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("HP301_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("HP301_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("HP301_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("HP301_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("HP301_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("HP301_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("HP301_Tzi8_all_pair.txt", sep = "\t")
Filtered_HP301_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_HP301_pan, file="Filtered_HP301_pan.csv")
# prepare files for the 12th genome Il14H 
# joint dataframe from previous two steps, watch out #variable change 
pan11 <- rbind(pan10,Filtered_HP301_pan[,-1])
post_pan11_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan11$Il14H), ";"))))
Il14H_pan_id <- as.matrix(Il14Hc[ ! Il14Hc$Query_gene %in% post_pan11_remove_list[,1],] )
colnames(Il14H_pan_id) <- "Query_gene"
write.csv(Il14H_pan_id, file = "Il14H_pan_id.txt")


# adding the 12th NAM genome: Il14H
# pan 12

Pan_reference <- read.csv("Il14H_pan_id.txt")
Pan_reference[,2]

Il14H <- Il14H_self
B73 <- read.csv("Il14H_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Il14H_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Il14H_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Il14H_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Il14H_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Il14H_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Il14H_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Il14H_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Il14H_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Il14H_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Il14H_HP301_all_pair.txt", sep = "\t")
P39 <- read.csv("Il14H_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Il14H_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Il14H_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Il14H_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Il14H_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Il14H_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Il14H_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Il14H_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Il14H_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Il14H_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Il14H_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Il14H_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Il14H_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Il14H_Tzi8_all_pair.txt", sep = "\t")
Filtered_Il14H_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Il14H_pan, file="Filtered_Il14H_pan.csv")
# prepare files for the 13th genome P39 
# joint dataframe from previous two steps, watch out #variable change
pan12 <- rbind(pan11,Filtered_Il14H_pan[,-1])
post_pan12_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan12$P39), ";"))))
P39_pan_id <- as.matrix(P39c[ ! P39c$Query_gene %in% post_pan12_remove_list[,1],] )
colnames(P39_pan_id) <- "Query_gene"
write.csv(P39_pan_id, file = "P39_pan_id.txt")

# adding the 13th NAM genome: P39
# pan 13

Pan_reference <- read.csv("P39_pan_id.txt")
Pan_reference[,2]

P39 <- P39_self
B73 <- read.csv("P39_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("P39_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("P39_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("P39_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("P39_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("P39_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("P39_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("P39_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("P39_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("P39_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("P39_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("P39_Il14H_all_pair.txt", sep = "\t")
CML52 <- read.csv("P39_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("P39_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("P39_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("P39_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("P39_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("P39_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("P39_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("P39_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("P39_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("P39_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("P39_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("P39_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("P39_Tzi8_all_pair.txt", sep = "\t")
Filtered_P39_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_P39_pan, file="Filtered_P39_pan.csv")
# prepare files for the 14th genome CML52 
# joint dataframe from previous two steps, watch out #variable change
pan13 <- rbind(pan12,Filtered_P39_pan[,-1])
post_pan13_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan13$CML52), ";"))))
CML52_pan_id <- as.matrix(CML52c[ ! CML52c$Query_gene %in% post_pan13_remove_list[,1],] )
colnames(CML52_pan_id) <- "Query_gene"
write.csv(CML52_pan_id, file = "CML52_pan_id.txt")


# adding the 14th NAM genome: CML52
# pan 14
#CML52
Pan_reference <- read.csv("CML52_pan_id.txt")
Pan_reference[,2]

CML52 <- CML52_self
B73 <- read.csv("CML52_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("CML52_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("CML52_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("CML52_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("CML52_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("CML52_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("CML52_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("CML52_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("CML52_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("CML52_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("CML52_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("CML52_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("CML52_P39_all_pair.txt", sep = "\t")
CML69 <- read.csv("CML52_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("CML52_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("CML52_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("CML52_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("CML52_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("CML52_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("CML52_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("CML52_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("CML52_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("CML52_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("CML52_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("CML52_Tzi8_all_pair.txt", sep = "\t")
Filtered_CML52_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_CML52_pan, file="Filtered_CML52_pan.csv")
# prepare files for the 15th genome CML69 
# joint dataframe from previous two steps, watch out #variable change
pan14 <- rbind(pan13,Filtered_CML52_pan[,-1])
post_pan14_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan14$CML69), ";"))))
CML69_pan_id <- as.matrix(CML69c[ ! CML69c$Query_gene %in% post_pan14_remove_list[,1],] )
colnames(CML69_pan_id) <- "Query_gene"
write.csv(CML69_pan_id, file = "CML69_pan_id.txt")


# adding the 15th NAM genome: CML69
# pan 15
#CML69

Pan_reference <- read.csv("CML69_pan_id.txt")
Pan_reference[,2]
CML69 <- CML69_self
B73 <- read.csv("CML69_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("CML69_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("CML69_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("CML69_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("CML69_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("CML69_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("CML69_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("CML69_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("CML69_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("CML69_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("CML69_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("CML69_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("CML69_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("CML69_CML52_all_pair.txt", sep = "\t")
CML103 <- read.csv("CML69_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("CML69_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("CML69_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("CML69_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("CML69_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("CML69_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("CML69_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("CML69_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("CML69_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("CML69_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("CML69_Tzi8_all_pair.txt", sep = "\t")
Filtered_CML69_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_CML69_pan, file="Filtered_CML69_pan.csv")
# prepare files for the 16th genome Ki11  
# joint dataframe from previous two steps, watch out #variable change
pan15 <- rbind(pan14,Filtered_CML69_pan[,-1])
post_pan15_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan15$Ki11), ";"))))
Ki11_pan_id <- as.matrix(Ki11c[ ! Ki11c$Query_gene %in% post_pan15_remove_list[,1],] )
colnames(Ki11_pan_id) <- "Query_gene"
write.csv(Ki11_pan_id, file = "Ki11_pan_id.txt")


# adding the 16th NAM genome: Ki11
# pan 16

Pan_reference <- read.csv("Ki11_pan_id.txt")
Pan_reference[,2]
Ki11 <- Ki11_self
B73 <- read.csv("Ki11_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Ki11_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Ki11_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Ki11_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Ki11_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Ki11_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Ki11_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Ki11_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Ki11_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Ki11_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Ki11_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Ki11_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Ki11_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Ki11_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Ki11_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Ki11_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Ki11_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Ki11_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Ki11_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Ki11_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Ki11_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Ki11_Ki3_all_pair.txt", sep = "\t")
NC350 <- read.csv("Ki11_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Ki11_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Ki11_Tzi8_all_pair.txt", sep = "\t")
Filtered_Ki11_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Ki11_pan, file="Filtered_Ki11_pan.csv")
# prepare files for the 17th genome CML228  
# joint dataframe from previous two steps, watch out #variable change
pan16 <- rbind(pan15,Filtered_Ki11_pan[,-1])
post_pan16_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan16$CML228), ";"))))
CML228_pan_id <- as.matrix(CML228c[ ! CML228c$Query_gene %in% post_pan16_remove_list[,1],] )
colnames(CML228_pan_id) <- "Query_gene"
write.csv(CML228_pan_id, file = "CML228_pan_id.txt")


# adding the 17th NAM genome: CML228
# pan 17
#CML228 
Pan_reference <- read.csv("CML228_pan_id.txt")
Pan_reference[,2]
CML228 <- CML228_self
B73 <- read.csv("CML228_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("CML228_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("CML228_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("CML228_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("CML228_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("CML228_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("CML228_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("CML228_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("CML228_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("CML228_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("CML228_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("CML228_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("CML228_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("CML228_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("CML228_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("CML228_CML103_all_pair.txt", sep = "\t")
CML247 <- read.csv("CML228_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("CML228_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("CML228_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("CML228_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("CML228_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("CML228_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("CML228_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("CML228_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("CML228_Tzi8_all_pair.txt", sep = "\t")
Filtered_CML228_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_CML228_pan, file="Filtered_CML228_pan.csv")
# prepare files for the 18th genome CML247  
# joint dataframe from previous two steps, watch out #variable change
pan17 <- rbind(pan16,Filtered_CML228_pan[,-1])
post_pan17_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan17$CML247), ";"))))
CML247_pan_id <- as.matrix(CML247c[ ! CML247c$Query_gene %in% post_pan17_remove_list[,1],] )
colnames(CML247_pan_id) <- "Query_gene"
write.csv(CML247_pan_id, file = "CML247_pan_id.txt")


# adding the 18th NAM genome: CML247
# pan 18

Pan_reference <- read.csv("CML247_pan_id.txt")
Pan_reference[,2]
CML247 <- CML247_self
B73 <- read.csv("CML247_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("CML247_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("CML247_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("CML247_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("CML247_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("CML247_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("CML247_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("CML247_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("CML247_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("CML247_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("CML247_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("CML247_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("CML247_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("CML247_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("CML247_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("CML247_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("CML247_CML228_all_pair.txt", sep = "\t")
CML277 <- read.csv("CML247_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("CML247_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("CML247_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("CML247_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("CML247_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("CML247_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("CML247_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("CML247_Tzi8_all_pair.txt", sep = "\t")
Filtered_CML247_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_CML247_pan, file="Filtered_CML247_pan.csv")
# prepare files for the 19th genome CML277   
# joint dataframe from previous two steps, watch out #variable change
pan18 <- rbind(pan17,Filtered_CML247_pan[,-1])
post_pan18_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan18$CML277), ";"))))
CML277_pan_id <- as.matrix(CML277c[ ! CML277c$Query_gene %in% post_pan18_remove_list[,1],] )
colnames(CML277_pan_id) <- "Query_gene"
write.csv(CML277_pan_id, file = "CML277_pan_id.txt")

# adding the 19th NAM genome: CML277
# pan 19
Pan_reference <- read.csv("CML277_pan_id.txt")
Pan_reference[,2]

CML277 <- CML277_self
B73 <- read.csv("CML277_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("CML277_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("CML277_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("CML277_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("CML277_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("CML277_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("CML277_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("CML277_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("CML277_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("CML277_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("CML277_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("CML277_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("CML277_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("CML277_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("CML277_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("CML277_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("CML277_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("CML277_CML247_all_pair.txt", sep = "\t")
CML322 <- read.csv("CML277_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("CML277_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("CML277_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("CML277_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("CML277_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("CML277_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("CML277_Tzi8_all_pair.txt", sep = "\t")
Filtered_CML277_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_CML277_pan, file="Filtered_CML277_pan.csv")
# prepare files for the 20th genome CML322
# joint dataframe from previous two steps, watch out #variable change
pan19 <- rbind(pan18,Filtered_CML277_pan[,-1])
post_pan19_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan19$CML322), ";"))))
CML322_pan_id <- as.matrix(CML322c[ ! CML322c$Query_gene %in% post_pan19_remove_list[,1],] )
colnames(CML322_pan_id) <- "Query_gene"
write.csv(CML322_pan_id, file = "CML322_pan_id.txt")



# adding the 20th NAM genome: CML322
# pan 20
# CML322
Pan_reference <- read.csv("CML322_pan_id.txt")
Pan_reference[,2]
CML322 <- CML322_self
B73 <- read.csv("CML322_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("CML322_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("CML322_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("CML322_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("CML322_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("CML322_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("CML322_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("CML322_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("CML322_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("CML322_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("CML322_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("CML322_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("CML322_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("CML322_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("CML322_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("CML322_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("CML322_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("CML322_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("CML322_CML277_all_pair.txt", sep = "\t")
CML333 <- read.csv("CML322_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("CML322_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("CML322_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("CML322_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("CML322_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("CML322_Tzi8_all_pair.txt", sep = "\t")
Filtered_CML322_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_CML322_pan, file="Filtered_CML322_pan.csv")
# prepare files for the 21th genome CML333
# joint dataframe from previous two steps, watch out #variable change
pan20 <- rbind(pan19,Filtered_CML322_pan[,-1])

write.csv(pan20, "~/Desktop/20_genomes_pan_rerun.csv") # this is the output for 20 genomes this time 

################don't run scripts below this point for now ###########################################





post_pan20_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan20$CML333), ";"))))
CML333_pan_id <- as.matrix(CML333c[ ! CML333c$Query_gene %in% post_pan20_remove_list[,1],] )
colnames(CML333_pan_id) <- "Query_gene"
write.csv(CML333_pan_id, file = "CML333_pan_id.txt")


# adding the 21th NAM genome: CML333
# pan 21
# CML333 

Pan_reference <- read.csv("CML333_pan_id.txt")
Pan_reference[,2]
CML333 <- CML333_self
B73 <- read.csv("CML333_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("CML333_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("CML333_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("CML333_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("CML333_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("CML333_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("CML333_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("CML333_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("CML333_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("CML333_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("CML333_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("CML333_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("CML333_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("CML333_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("CML333_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("CML333_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("CML333_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("CML333_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("CML333_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("CML333_CML322_all_pair.txt", sep = "\t")
Ki3 <- read.csv("CML333_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("CML333_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("CML333_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("CML333_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("CML333_Tzi8_all_pair.txt", sep = "\t")
Filtered_CML333_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_CML333_pan, file="Filtered_CML333_pan.csv")
# prepare files for the 22th genome Ki3
# joint dataframe from previous two steps, watch out #variable change
pan21 <- rbind(pan20,Filtered_CML333_pan[,-1])
post_pan21_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan21$Ki3), ";"))))
Ki3_pan_id <- as.matrix(Ki3c[ ! Ki3c$Query_gene %in% post_pan21_remove_list[,1],] )
colnames(Ki3_pan_id) <- "Query_gene"
write.csv(Ki3_pan_id, file = "Ki3_pan_id.txt")


# adding the 22th NAM genome: Ki3
# pan 22
# Ki3 

Pan_reference <- read.csv("Ki3_pan_id.txt")
Pan_reference[,2]

Ki3 <- Ki3_self
B73 <- read.csv("Ki3_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Ki3_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Ki3_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Ki3_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Ki3_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Ki3_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Ki3_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Ki3_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Ki3_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("Ki3_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("Ki3_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Ki3_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Ki3_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Ki3_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Ki3_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Ki3_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Ki3_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Ki3_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Ki3_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Ki3_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Ki3_CML333_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Ki3_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Ki3_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Ki3_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Ki3_Tzi8_all_pair.txt", sep = "\t")
Filtered_Ki3_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Ki3_pan, file="Filtered_Ki3_pan.csv")
# prepare files for the 23th genome CML103
# joint dataframe from previous two steps, watch out #variable change
pan22 <- rbind(pan21,Filtered_Ki3_pan[,-1])
post_pan22_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan22$CML103), ";"))))
CML103_pan_id <- as.matrix(CML103c[ ! CML103c$Query_gene %in% post_pan22_remove_list[,1],] )
colnames(CML103_pan_id) <- "Query_gene"
write.csv(CML103_pan_id, file = "CML103_pan_id.txt")


# adding the 23th NAM genome: CML103
# pan 23
#CML103 
Pan_reference <- read.csv("CML103_pan_id.txt")
Pan_reference[,2]

CML103 <- CML103_self
B73 <- read.csv("CML103_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("CML103_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("CML103_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("CML103_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("CML103_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("CML103_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("CML103_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("CML103_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("CML103_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("CML103_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("CML103_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("CML103_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("CML103_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("CML103_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("CML103_CML69_all_pair.txt", sep = "\t")
CML228 <- read.csv("CML103_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("CML103_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("CML103_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("CML103_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("CML103_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("CML103_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("CML103_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("CML103_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("CML103_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("CML103_Tzi8_all_pair.txt", sep = "\t")
Filtered_CML103_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_CML103_pan, file="Filtered_CML103_pan.csv")
# prepare files for the 24th genome Tx303 
# joint dataframe from previous two steps, watch out #variable change
pan23 <- rbind(pan22,Filtered_CML103_pan[,-1])
post_pan23_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan23$Tx303), ";"))))
Tx303_pan_id <- as.matrix(Tx303c[ ! Tx303c$Query_gene %in% post_pan23_remove_list[,1],] )
colnames(Tx303_pan_id) <- "Query_gene"
write.csv(Tx303_pan_id, file = "Tx303_pan_id.txt")

# adding the 24th NAM genome: Tx303
# pan 24

Pan_reference <- read.csv("Tx303_pan_id.txt")
Pan_reference[,2]

Tx303 <- Tx303_self 
B73 <- read.csv("Tx303_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("Tx303_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("Tx303_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("Tx303_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("Tx303_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("Tx303_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("Tx303_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("Tx303_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("Tx303_Mo18W_all_pair.txt", sep = "\t")
HP301 <- read.csv("Tx303_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("Tx303_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("Tx303_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("Tx303_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("Tx303_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("Tx303_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("Tx303_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("Tx303_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("Tx303_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("Tx303_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("Tx303_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("Tx303_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("Tx303_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("Tx303_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("Tx303_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("Tx303_Tzi8_all_pair.txt", sep = "\t")
Filtered_Tx303_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_Tx303_pan, file="Filtered_Tx303_pan.csv")
# prepare files for the 25th genome NC358 
# joint dataframe from previous two steps, watch out #variable change
pan24 <- rbind(pan23,Filtered_Tx303_pan[,-1])
post_pan24_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan24$NC358), ";"))))
NC358_pan_id <- as.matrix(NC358c[ ! NC358c$Query_gene %in% post_pan24_remove_list[,1],] )
colnames(NC358_pan_id) <- "Query_gene"
write.csv(NC358_pan_id, file = "NC358_pan_id.txt")

# adding the 25th NAM genome: NC358
# pan 25

Pan_reference <- read.csv("NC358_pan_id.txt")
Pan_reference[,2]

NC358 <- NC358_self
B73 <- read.csv("NC358_B73_all_pair.txt", sep = "\t")
B97 <- read.csv("NC358_B97_all_pair.txt", sep = "\t")
Ky21 <- read.csv("NC358_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("NC358_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("NC358_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("NC358_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("NC358_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("NC358_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("NC358_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("NC358_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("NC358_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("NC358_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("NC358_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("NC358_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("NC358_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("NC358_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("NC358_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("NC358_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("NC358_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("NC358_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("NC358_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("NC358_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("NC358_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("NC358_NC350_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("NC358_Tzi8_all_pair.txt", sep = "\t")
Filtered_NC358_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_NC358_pan, file="Filtered_NC358_pan.csv")
# prepare files for the 26th genome B97
# joint dataframe from previous two steps, watch out #variable change
pan25 <- rbind(pan24,Filtered_NC358_pan[,-1])
post_pan25_remove_list <- as.matrix(unique(unlist(strsplit(as.character(pan25$B97), ";"))))
B97_pan_id <- as.matrix(B97c[ ! B97c$Query_gene %in% post_pan25_remove_list[,1],] )
colnames(B97_pan_id) <- "Query_gene"
write.csv(B97_pan_id, file = "B97_pan_id.txt")


# adding the 26th NAM genome: B97
# pan 26

Pan_reference <- read.csv("B97_pan_id.txt")
Pan_reference[,2]
B97 <- B97_self
B73 <- read.csv("B97_B73_all_pair.txt", sep = "\t")
Ky21 <- read.csv("B97_Ky21_all_pair.txt", sep = "\t")
M162W <- read.csv("B97_M162W_all_pair.txt", sep = "\t")
Ms71 <- read.csv("B97_Ms71_all_pair.txt", sep = "\t")
Oh7B <- read.csv("B97_Oh7B_all_pair.txt", sep = "\t")
Oh43 <- read.csv("B97_Oh43_all_pair.txt", sep = "\t")
M37W <- read.csv("B97_M37W_all_pair.txt", sep = "\t")
Mo18W <- read.csv("B97_Mo18W_all_pair.txt", sep = "\t")
Tx303 <- read.csv("B97_Tx303_all_pair.txt", sep = "\t")
HP301 <- read.csv("B97_HP301_all_pair.txt", sep = "\t")
Il14H <- read.csv("B97_Il14H_all_pair.txt", sep = "\t")
P39 <- read.csv("B97_P39_all_pair.txt", sep = "\t")
CML52 <- read.csv("B97_CML52_all_pair.txt", sep = "\t")
CML69 <- read.csv("B97_CML69_all_pair.txt", sep = "\t")
CML103 <- read.csv("B97_CML103_all_pair.txt", sep = "\t")
CML228 <- read.csv("B97_CML228_all_pair.txt", sep = "\t")
CML247 <- read.csv("B97_CML247_all_pair.txt", sep = "\t")
CML277 <- read.csv("B97_CML277_all_pair.txt", sep = "\t")
CML322 <- read.csv("B97_CML322_all_pair.txt", sep = "\t")
CML333 <- read.csv("B97_CML333_all_pair.txt", sep = "\t")
Ki3 <- read.csv("B97_Ki3_all_pair.txt", sep = "\t")
Ki11 <- read.csv("B97_Ki11_all_pair.txt", sep = "\t")
NC350 <- read.csv("B97_NC350_all_pair.txt", sep = "\t")
NC358 <- read.csv("B97_NC358_all_pair.txt", sep = "\t")
Tzi8 <- read.csv("B97_Tzi8_all_pair.txt", sep = "\t")
Filtered_B97_pan <- left_join(Pan_reference, B73) %>% left_join(Tzi8) %>% left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>% left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97)
write.csv(Filtered_B97_pan, file="Filtered_B97_pan.csv")

#merging everything together, make Pan 26
pan26 <- rbind(pan25,Filtered_B97_pan[-1])
write.csv(pan26, file ="~/Desktop/NAM_26_genomes/nucmer_1000_10_sub_genome/pan26_all.csv")

