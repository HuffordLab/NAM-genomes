library(abc)
library(tidyverse)
par(ask=F)

#set.seed(123)
con <- file("~/Dropbox/NAM_DFE_ABC/pop.txt","r")
first_line <- readLines(con,n=1)
close(con)
col_string <- str_split(first_line, " ")[[1]] %>% 
  str_remove(pattern = "#") %>% 
  stringi::stri_remove_empty()

full_df <- read_delim(
  "~/Dropbox/NAM_DFE_ABC/pop.txt", delim = " ",
  col_names = col_string, comment = "#") %>% 
  drop_na() 

  #some TD values are NA
  # select(-c("bin_mean1",      "bin_mean_dist1", "wt_mean_dist1",  "bin_sd_dist1",   "wt_sd_dist1" ,
  #           "bin_mean2",      "bin_mean_dist2", "wt_mean_dist2",  "bin_sd_dist2",   "wt_sd_dist2" ,
  #           "bin_mean3",      "bin_mean_dist3", "wt_mean_dist3",  "bin_sd_dist3",   "wt_sd_dist3",
  #           "bin_mean4",      "bin_mean_dist4", "wt_mean_dist4",  "bin_sd_dist4",   "wt_sd_dist4"
  # 
  #           ))


# select(-c(
#   "sfs_pi1", "sfs_td1","sfs_pi2", "sfs_td2", "sfs_pi3","sfs_td3"
# ))
names(full_df)

mean(full_df$Na > full_df$Nb)
plot(full_df$Na, full_df$Nb)

pc_sumstat_df <- prcomp(
  drop_na(select(full_df, 
         -c(Na, N0, Nb, B_t, 
            sfs1_shape, sfs1_mean, 
            sfs2_shape, sfs2_mean, 
            p_neutral, p_sfs1)
  )),
  center = F, scale. = T)
plot(pc_sumstat_df$x[,1], pc_sumstat_df$x[,2])
(pc_sumstat_df$sdev^2/sum(pc_sumstat_df$sdev^2))

varz <- cumsum((pc_sumstat_df$sdev ^ 2) / sum(pc_sumstat_df$sdev ^ 2)) <= 1
mx_col <- max(c(2,length(varz[varz])))
pc_sumstat_df <- data.frame(pc_sumstat_df$x) %>% select(1:mx_col)


hist(log(full_df$sfs_pi1))

abline(v = log(0.008))
hist(full_df$sfs_pi1)
abline(v = mean(full_df$sfs_pi1), col = "red")
abline(v = 0.008)
summary(full_df$sfs_pi1)
hist(full_df$sfs_td1)
range(full_df$sfs_td1)
round(quantile(full_df$sfs_pi1, c(0.025, 0.975)), 3)
round(quantile(full_df$sfs_td1, c(0.025, 0.975)), 3)
#pairs(select(full_df, -c(Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1)))
corrs <- cor(full_df)
plot(density(corrs))

target_idx <- sample(seq_len(nrow(full_df)), 1)
target_df <- full_df[target_idx, ]
full_df <- full_df[-target_idx, ]

sumstat_df <- select(full_df, -c(Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1))
param_df <- select(full_df, c(Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1))

pc_df <- data.frame(pc_sumstat_df) %>% slice(-target_idx)
pc_target <- data.frame(pc_sumstat_df) %>% slice(target_idx)


target_stats <- select(target_df, -c(Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1))
target_params <- select(target_df, c(Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1))


res.gfit.bott <- gfit(target=target_stats, 
                   sumstat=sumstat_df,
                   statistic=mean, nb.replicate=10)
summary(res.gfit.bott)$pvalue
 
# res <- abc(target=target_stats,
#            param=param_df,
#            sumstat=sumstat_df,
#            tol=0.01, transf=c("log"), method="neuralnet")

res <- abc(target=pc_target,
           param=param_df,
           sumstat=pc_df,
           tol=0.01, transf=c("log"), method="neuralnet")

posts_df <- data.frame(res$adj.values)
names(posts_df)
c_name <- "sfs1_shape"
plot(density(posts_df[[c_name]]))

posts_df[[c_name]]

den_prior <- density(param_df[[c_name]])
den_post <- density(posts_df[[c_name]])
plot(den_prior, 
     xlim = range(c(range(den_prior$x), range(den_post$x))),
     ylim = range(c(range(den_prior$y), range(den_post$y))),
     col = "red"
)
abline(v = mean(param_df[[c_name]]), col = "red")

lines(den_post, col = "blue")
abline(v = target_df[[c_name]])
abline(v = mean(posts_df[[c_name]]), col = "blue")
abline(v = target_params[[c_name]])
diff(quantile(param_df[[c_name]], c(0.025, 0.975))) > diff(quantile(posts_df[[c_name]], c(0.025, 0.975)))

quantile(param_df[[c_name]], c(0.025, 0.975))
quantile(posts_df[[c_name]], c(0.025, 0.975))

mean(posts_df[[c_name]] > target_params[[c_name]])
posts_df <- data.frame(res$unadj.values)
mean(posts_df[[c_name]] > target_params[[c_name]])

