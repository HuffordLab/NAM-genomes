library(abc)
library(tidyverse)
library(optparse)
par(ask=F)

option_list <-  list(
  make_option(c("-s", "--sim_sfs"), type="character", default=NULL, 
              help="Simulated SFS file. Space delimited", metavar="character"),
  make_option(c("-p", "--sim_param"), type="character", default=NULL, 
              help="Simulated file listing true parameters. Space delimited", metavar="character"),
  make_option(c("-e", "--obs_sfs"), type="character", default=NULL, 
              help="Observed sfs file, each line contains SFS counts for 4fold, 0fold, and sv categories.", metavar="character"),
  make_option(c("-c", "--sfs_count"), type="numeric", default=0, 
              help="Minimum raw count of required for window to be analyzed", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Name of output file containing all posterior draws.", metavar="character")
)

opt_parser <-  OptionParser(option_list=option_list)
opt <-  parse_args(opt_parser)

n_ind <- 26
fold_bins <- (1:(n_ind/2))
to_keep <- c(fold_bins, fold_bins + n_ind*1, fold_bins + n_ind*2)

#-----------------------
#simulated SFS
full_df <- vroom::vroom(opt$sim_sfs, delim = " ", col_names = FALSE) %>%
  dplyr::select(-ncol(.)) %>% 
  dplyr::select(all_of(to_keep)) %>%
  mutate(sm  = rowSums(.)) %>%
  ungroup() %>% 
  rowwise() %>% 
  mutate_all(.funs = ~ .x / sm) %>% 
  dplyr::select(-sm)

#3 observed SFS types pasted together
#sv_sum <- apply(read_delim(opt$obs_sfs, delim = "\t", col_names = FALSE)[ fold_bins + n_ind*2], 1, sum) > opt$sfs_count
sv_sum <- apply(read_delim(opt$obs_sfs, delim = "\t", col_names = FALSE), 1, sum) > opt$sfs_count
window_targets <- which(sv_sum)
print("window targets")
print(window_targets)
#compute ALL windows
#window_tartets <- seq_along(sv_sum)

if (sum(sv_sum) > 0){ #filter for empty data set
  #observed nam sfs data
  nam_df <- vroom::vroom(
    opt$obs_sfs, delim = "\t",
    col_names = FALSE
  ) %>%
    dplyr::select(all_of(to_keep)) %>% 
    mutate(sm  = rowSums(.)) %>%
    rowwise() %>% 
    mutate_all(.funs = ~ .x / sm) %>% 
    dplyr::select(-sm) %>% 
    ungroup()


 params <- 
  c("Na", "N0", "Nb", "B_t", "mu_sv",
    "sfs1_shape", "sfs1_mean", "sfs2_shape", "sfs2_mean")
  
  
  #-----------------------
  #simulated parameters
  param_df <- read_delim(
    opt$sim_param,
    col_names = params,
    delim = " "
  ) %>% 
    dplyr::select(-c(Na, B_t, mu_sv))
  
  
  all_abc_post <- 
    window_targets %>%  #!!!!
    #1:nrow(nam_df) %>%
    map_df(function(target_idx){
      #raw
      target_stats <- nam_df[target_idx, ]
      res <- abc(target=target_stats,
                 param=param_df,
                 sumstat=full_df,
                 tol=0.002,
                 transf=c("none"),
                 method = "neuralnet",
                 sizenet = 2)
      
      posts_df <- data.frame(res$adj.values) 
      print(posts_df[1,])
      #mean_post <- summarise_all(posts_df, .funs = median)
      #prior_low <- summarise_all(param_df, .funs = min)
      #prior_high <- summarise_all(param_df, .funs = max)
      #post_in <- mean(mean_post > prior_low & mean_post < prior_high)
      
      #if (post_in != "1") {
      fname = str_split(str_split(opt$obs_sfs, "NAM_SFS[_-]")[[1]][2], "_all")[[1]][1]
      posts_df %>% mutate(window = target_idx, sfs = fname)
      #} else {
      #  tibble()
      #}
    })
  
  #-----------------------
  write_csv(all_abc_post, opt$out)
} else {
  print("No windows have sufficient SV data. Writing empty tibble to file")
  write_csv(tibble(), opt$out)
}


