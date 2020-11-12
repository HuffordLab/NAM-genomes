library(abc)
library(tidyverse)
library(ggridges)
library(patchwork)
library(cowplot)
library(magrittr)
library(purrr)
theme_set(theme_cowplot())
par(ask=F)

n_ind <- 26
fold_bins <- (1:(n_ind/2))
to_keep <- c(fold_bins, fold_bins + n_ind*1, fold_bins + n_ind*2)
#var0 <- function(x) var(x) > 0


#posterior simulated
post_sfs <- read_delim(
  "predict/data/postcheck_out/post_sfs.txt",
  delim = " ",
  col_names = FALSE
) %>%
  dplyr::select(-ncol(.)) %>% 
  dplyr::select(all_of(to_keep)) 

#empirical
emprical_sfs <- read_delim(
  "predict/data/sfs/NAM_SFS-dup_inv_tra_ins_del_all.txt", delim = "\t",
  col_names = FALSE
) %>%
  dplyr::select(-ncol(.)) %>% 
  dplyr::select(all_of(to_keep)) #%>% 
# mutate(sm  = rowSums(.)) %>%
# rowwise() %>%
# mutate_all(.funs = ~ .x / sm) %>%
# dplyr::select(-sm) %>% 
# ungroup()

post_params <- read_delim(
  "predict/data/postcheck_out/post_params.txt",
  col_names = TRUE,
  delim = "\t"
) %>% 
  mutate(idx = 1:n())

posts_df <- 1:nrow(emprical_sfs) %>% 
  map_df(function(wind){
    
    #turn 3 sfs for nam data for current window into a df
    emp_wind <- tibble(empirical = as.vector(t(slice(emprical_sfs, wind)))) %>% 
      mutate(idx = 1:n())
    
    #match rows of simulated sfs to parmeter df to get rows for the current window 
    post_draw_wind <- unname(t(slice(post_sfs, pull(filter(post_params, window == wind), idx)))) 
    
    #in case no data is present for a window
    if(length(post_draw_wind) < 1){
      tibble()
    } else {
      
    #make workable names and add index
    post_draw_wind <- set_colnames(x = post_draw_wind, paste0("x", 1:ncol(post_draw_wind))) %>% 
      as_tibble() %>% 
      mutate(idx = 1:n())
    #join the emprical and simulated data, add all simulated data to single row
    post_df <- full_join(emp_wind, post_draw_wind, by = "idx") %>% 
      mutate(type = rep(c("4-fold", "0-fold", "SV"), each = 13)) %>% 
      pivot_longer(cols = starts_with("x"), names_to = "draw", values_to = "simulated") %>% 
      mutate(window = wind,
             mu_sv = post_params$musv[wind])
    }
  })

posts_df$type <- factor(posts_df$type, levels = c("4-fold", "0-fold", "SV"))

a <- posts_df %>% 
  ggplot(aes(empirical, simulated, group = draw)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  facet_wrap(~type, scales = "free") 
  #facet_wrap(~type + mu_sv, scales = "free") 
a


max(posts_df$simulated - posts_df$empirical)
max(posts_df$empirical)
max(posts_df$simulated)

b <- posts_df %>%
  #filter(draw == sample(posts_df$draw, 1),
  #       window ==sample(posts_df$window, 1)) %>% 
  ggplot(aes(idx, (simulated - empirical)/(empirical+1), group = window)) +
  geom_jitter(height = 0, width = 0.25, alpha = 0.05) +
  geom_smooth(aes(idx, (simulated - empirical)/(empirical+1)), inherit.aes = F, colour = "dodgerblue") +
  geom_hline(yintercept = 0, lty = 2, colour = "forestgreen", lwd = 1.2) +
  facet_wrap(~type, scales = "free")  
  #facet_wrap(~type + mu_sv, scales = "free")  
b

b <- posts_df %>%
  #filter(draw == sample(posts_df$draw, 1),
  #       #window ==sample(posts_df$window, 1)
  #       ) %>% 
  ggplot(aes(idx, (simulated - empirical)/(empirical+1), group = idx)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_smooth(aes(idx, (simulated - empirical)/(empirical+1)), inherit.aes = F, colour = "dodgerblue") +
  #geom_hline(yintercept = 0, lty = 2, colour = "forestgreen", lwd = 1.2) +
  ylim(0, 200) +
  facet_wrap(~type, scales = "free_x")


c <- posts_df %>% 
  ggplot(aes(idx, empirical, group = window)) +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(idx, empirical), inherit.aes = F, colour = "dodgerblue") +
  facet_wrap(~type, nrow = 3, scales = "free") 
c

c + ggtitle("A") + xlab("allele frequency") + b + ggtitle("B") + xlab("allele frequency") + plot_layout(ncol = 2, widths = c(0.2, 0.8)) + 
  ggsave(filename = "analysis_viz/figures/postcheck_viz.pdf", width = 15, height = 10)


posts_df %>% 
  group_by(type, mu_sv) %>% 
  summarise(prop_gt = mean(empirical >= simulated)) %>% 
  arrange(mu_sv, type)
