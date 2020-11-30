#Permutation Test 
#Are the expression candidates enriched for significant GWAS SVs?

library(tidyverse)
library(ggplot2)

#GWAS hits: `NAM_SV_T1.mlma`
#candidate SV coordinates: `all_candSV_coordinates.bed`
#candidate SV formated to naming of GWAS file: `list.txt`
#overlap between list and GWAS hits: `SV_GWAS_hit_overlap.mlma`

candSV_GWAS_hits<-read_tsv("SV_GWAS_hit_overlap.mlma", col_names = c("chromosome","SNP","physical_position","ref_allele","alt_allele","ref_allele_freq","SNP_effect","std_error","p_value"))
Full_GWAS_hits<-read_tsv("NAM_SV_T1.mlma")

filter(candSV_GWAS_hits, p_value <= 0.05) %>% nrow()

run_perm<-function(input, l){ #input = reference p-value column, l=sample size you want to draw i.e.length of candidate list
  sim_sig_p<-vector() #empty vector for results
  counter=0 #Sets the counter at 0
  while(counter < 1000){ #starts the while loop which let's it iterate 1000 times
    run<-sample(input,l)  #Takes a random sample from reference the size of l
    num_sig<-run %>% .[. <= 0.05] %>% length() #calculates the number of significant hits in the sample
    sim_sig_p<-c(sim_sig_p,num_sig) #adds results to vector
    counter=(counter+1) #moves the counter forward by 1
  }
  nm<-tibble(sim_sig_p)
  return(nm)
}

perm_results<-run_perm(Full_GWAS_hits$p, nrow(candSV_GWAS_hits))
ggplot(perm_results, aes(x=sim_sig_p)) +
  geom_histogram(binwidth = 1, stat = "count") +
  geom_vline(xintercept = 3, color= "red") #number of our candidates that are significant

filter(perm_results, sim_sig_p < 3) %>% nrow() #number of permutations with fewer significant results
filter(perm_results, sim_sig_p > 3) %>% nrow() #number of permutations with more significant results
filter(perm_results, sim_sig_p == 3) %>% nrow() #number of permutations with the same number of significant results

#Trying the permutation test but with median p-values
run_perm_med<-function(input, l){ #input = reference p-value column, l=sample size you want to draw i.e.length of candidate list
  sim_med_p<-vector() #empty vector for results
  counter=0 #Sets the counter at 0
  while(counter < 1000){ #starts the while loop which let's it iterate 1000 times
    run<-sample(input,l)  #Takes a random sample from reference the size of l
    p_med<-run %>% median(na.rm = TRUE) #calculates the median of the p-value
    sim_med_p<-c(sim_med_p,p_med) #adds results to vector
    counter=(counter+1) #moves the counter forward by 1
  }
  nm<-tibble(sim_med_p)
  return(nm)
}

perm_results_med<-run_perm_med(Full_GWAS_hits$p, nrow(candSV_GWAS_hits))
candSV_GWAS_hits$p_value %>% median(na.rm = TRUE)
ggplot(perm_results_med, aes(x=sim_med_p)) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept = 0.406743, color= "red") #number of our candidates that are significant

filter(perm_results_med, sim_med_p < 0.406743) %>% nrow() #number of permutations with lower median p-values
filter(perm_results_med, sim_med_p > 0.406743) %>% nrow() #number of permutations with higher median p-values
filter(perm_results_med, sim_med_p == 0.406743) %>% nrow() #number of permutations with the same median p-value
