library(abc)
library(tidyverse)
library(ggridges)
library(patchwork)
library(cowplot)
library(magrittr)
theme_set(theme_cowplot(font_size = 32))
par(ask=F)

n_ind <- 26
fold_bins <- (1:(n_ind/2))
to_keep <- c(fold_bins, fold_bins + n_ind*1, fold_bins + n_ind*2)
#var0 <- function(x) var(x) > 0


full_df <- vroom::vroom(
  "predict/pop_sfs.txt", delim = " ",
  #sfs_file, delim = " ",
  col_names = FALSE
) %>%
  dplyr::select(-ncol(.)) %>% 
  dplyr::select(all_of(to_keep)) %>% 
  mutate(sm  = rowSums(.)) %>%
  rowwise() %>%
  mutate_all(.funs = ~ .x / sm) %>%
  dplyr::select(-sm) %>% 
  ungroup()


params <- 
  c("Na", "N0", "Nb", "B_t", "mu_sv",
    "sfs1_shape", "sfs1_mean", "sfs2_shape", "sfs2_mean")

param_df_full <- read_delim(
  "predict/pop.txt",
  #params_file,
  col_names = params,
  delim = " "
) %>% 
  dplyr::select(-c(Na, B_t))


abc_val <- function(full_df, param_df_full){

  itts <- 200
  p <- progress_estimated(itts)
  set.seed(1133)
  
  val_df <- 
    sample(seq_len(nrow(full_df)), itts) %>% 
    map_df(function(target_idx){
      print(p$tick())
     
      #raw
      target_stats <- full_df[target_idx, ]
      sumstat_df <- full_df[-target_idx, ]
      
      target_df <- param_df_full[target_idx, ]
      param_df <-  param_df_full[-target_idx, ]
      
      
      res <- abc(target=target_stats,
                 param=param_df,
                 sumstat=sumstat_df,
                 tol=0.002,
                 transf=c("none"),
                 method = "neuralnet", 
                 sizenet = 2
      )
      
      posts_df <- data.frame(res$adj.values)
      
      names(posts_df) %>% 
        map_df(function(c_name){
          prop_gt <- mean(posts_df[[c_name]] >= target_df[[c_name]])
          creds <- quantile(posts_df[[c_name]], c(0,0.025,0.975,1))
          w.in_cred <- as.numeric(creds[2] < target_df[[c_name]] & creds[3] > target_df[[c_name]])
          w.in_prior <- as.numeric(creds[1] < target_df[[c_name]] & creds[4] > target_df[[c_name]])
          diff_sc <- mean( (posts_df[[c_name]] - target_df[[c_name]]) / (target_df[[c_name]] + 1) )
          var_sc <- sd(posts_df[[c_name]] / sd(param_df[[c_name]]))
          mean_diff <- mean(sqrt((posts_df[[c_name]] -  target_df[[c_name]])^2) /( target_df[[c_name]] + 1))
          mean_post <- mean(posts_df[[c_name]])
          tibble(mean_diff, parameter = c_name, w.in_cred, w.in_prior, prop_gt, diff_sc, var_sc, mean_post)
        }) %>% 
        bind_cols(. , target_df)
      
    })

  val_df %<>% 
    mutate(is_zero = factor(if_else(sfs1_mean == 0 | sfs2_mean == 0, true = "zero", false = "not_zero")))

  val_df
}


viz_val <- function(val_df, bns = 30){
  
  val_df %<>%
    mutate(
      parameter = str_replace(parameter, "sfs2_shape", "shape_sv"),
      parameter = str_replace(parameter, "sfs1_shape", "shape_0"),
      parameter = str_replace(parameter, "sfs2_mean", "s_sv"),
      parameter = str_replace(parameter, "sfs1_mean", "s_0"),
    )
  pgt <- val_df %>% 
    ggplot(aes(x = prop_gt, y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.01) +
    geom_vline(xintercept = 0.5, lty = 2) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme(legend.position = "n") 
  
  cred <- val_df %>% 
    ggplot(aes(x = w.in_cred, y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = 2, scale = 0.5, rel_min_height = 0.01) +
    theme(axis.text.y=element_blank(),
          legend.position = "n") +
    scale_x_continuous(breaks = c(0,1)) +
    ylab("")
  
  wpr <- val_df %>% 
    ggplot(aes(x = w.in_prior, y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = 2, scale = 0.5, rel_min_height = 0.01) +
    theme(axis.text.y=element_blank(),
          legend.position = "n") +
    scale_x_continuous(breaks = c(0,1)) +
    ylab("")
  
  
  df_sc <- val_df %>% 
    ggplot(aes(x = diff_sc, y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, lty = 2) +
    theme(axis.text.y=element_blank(),
          legend.position = "n") +
    ylab("")
  
  v_sc <- 
    val_df %>%
    filter(sfs1_mean != 0) %>% 
    ggplot(aes(x = log(var_sc), y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.001) +
    geom_vline(xintercept = log(1), lty = 2, lwd =1.2) +
    theme(axis.text.y=element_blank(),
          legend.position = "n") +
    ylab("")
  
  pgt + cred + wpr + v_sc + plot_layout(nrow = 1) 
}


val_df <- abc_val(full_df, param_df_full)
write_csv(x = val_df, path = "analysis_viz/val_mu.csv")

val_df <- read_csv("analysis_viz/val_mu.csv")

bns <- 30
viz_val(val_df) + ggsave("analysis_viz/figures/abc_validate.pdf", width = 12, height = 10)

v_sc + xlim(-5, 1)

cred <- val_df %>% 
  ggplot(aes(x = w.in_cred, y = parameter, fill = parameter, colour =  parameter)) +
  geom_density_ridges(stat = "binline", bins = 2, scale = 0.5, rel_min_height = 0.01) +
  theme(legend.position = "n") +
  scale_x_continuous(breaks = c(0,1)) +
  ylab("")

v_sc <- 
  val_df %>%
  ggplot(aes(x = log(var_sc), y = parameter, fill = parameter, colour =  parameter)) +
  geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.001) +
  geom_vline(xintercept = log(1), lty = 2, lwd =1.2) +
  theme(axis.text.y=element_blank(),
        legend.position = "n") +
  ylab("")

cred + v_sc + xlim(-5, 1) + plot_layout(nrow = 1)

val_df %>% 
  group_by(parameter) %>% 
  summarise(vs = mean(log(var_sc)),
            in_cred = mean(w.in_cred == 1)
  )


val_df %>% 
  select(parameter, mean_post)

val_df %>% 
  filter(parameter %in% c("sfs1_mean", "sfs2_mean")) %>%
  mutate(var = mean_post) %>%
  ggplot(aes(parameter, -var)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.05) +
  facet_grid(~is_zero)


val_df %>% 
  filter(parameter %in% c("sfs1_mean", "sfs2_mean")) %>%
  mutate(var = mean_post) %>%
  ggplot(aes(mu_sv, -var, colour = is_zero)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~parameter) 

val_df %>% 
  filter(parameter %in% c("sfs1_mean")) %>%
  mutate(var = mean_post) %>%
  ggplot(aes(mu_sv, -var, colour = is_zero)) +
  geom_point() +
  geom_smooth()


val_df %>% 
  filter(parameter %in% c("sfs2_mean", "mu_sv")) %>%
  group_by(parameter, is_zero) %>% 
  mutate(idx = 1:n()) %>% 
  dplyr::select(idx, parameter, mean_post) %>% 
  pivot_wider(id_cols = c(idx, is_zero), names_from = "parameter", values_from = mean_post) %>% 
  ggplot(aes(mu_sv, -sfs2_mean, colour = is_zero)) +
  geom_point() +
  geom_smooth(method = "lm")

val_df %>% 
  filter(parameter %in% c("sfs1_mean", "mu_sv")) %>%
  group_by(parameter, is_zero) %>% 
  mutate(idx = 1:n()) %>% 
  dplyr::select(idx, parameter, mean_post) %>% 
  pivot_wider(id_cols = c(idx, is_zero), names_from = "parameter", values_from = mean_post) %>% 
  ggplot(aes(mu_sv, -sfs1_mean, colour = is_zero)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(legend.position = "n")


val_df %>% 
  filter(parameter %in% c("sfs2_mean"), sfs2_mean == 0) %>% 
  lm(-mean_post ~ N0 + Nb + mu_sv + sfs2_shape + sfs1_shape, data = .) %>% 
  summary()

val_df %>% 
  filter(parameter %in% c("sfs1_mean"), sfs2_mean == 0) %>% 
  lm(-mean_post ~ N0 + Nb + mu_sv + sfs2_shape + sfs1_shape, data = .) %>% 
  summary()


val_df %>% 
  filter(parameter %in% c("sfs2_mean"), sfs2_mean != 0) %>% 
  lm(-mean_post ~ N0 + Nb + mu_sv + sfs2_shape + sfs1_shape, data = .) %>% 
  summary()

val_df %>% 
  filter(parameter %in% c("sfs1_mean"), sfs2_mean != 0) %>% 
  lm(-mean_post ~ N0 + Nb + mu_sv + sfs2_shape + sfs1_shape, data = .) %>% 
  summary()

val_df %>% 
  filter(sfs1_mean == 0,
         parameter %in% c("sfs2_mean", "sfs1_mean")) %>% 
  select(mean_post)


#spot checking


target_idx <- sample(which(param_df_full$sfs1_mean == 0), 1)
target_idx <- sample(which(param_df_full$sfs1_mean != 0), 1)
param_df_full[target_idx,]
#raw
target_stats <- full_df[target_idx, ]
sumstat_df <- full_df[-target_idx, ]

target_df <- param_df_full[target_idx, ]
param_df <-  param_df_full[-target_idx, ]


res <- abc(target=target_stats,
           param=param_df,
           sumstat=sumstat_df,
           tol=0.002,
           transf=c("none"),
           method = "neuralnet", 
           sizenet = 2)
summary(res)


posts_df <- data.frame(res$adj.values)
c_name <- "sfs2_mean"
mean(posts_df$sfs1_mean)
hist(posts_df$sfs1_mean)
hist(posts_df$sfs2_mean)

tibble(
  value = c(-param_df[[c_name]], -posts_df[[c_name]], 0),
  distribution = c(
    rep("prior", nrow(param_df)), 
    rep("posterior", nrow(posts_df)), 
    "simulated")) %>% 
  ggplot(aes(x = value, colour = distribution)) +
  geom_density(fill=NA, lwd = 2) +
  geom_vline(xintercept = -target_df[[c_name]], lwd = 2) +
  scale_colour_manual(values = c("prior" = "grey60", "posterior" = "blue", "simulated" = "black")) +
  xlab("sfs2_mean") +
  xlim(-0.11, 0.01) +
  theme_cowplot(font_size = 35)

log(sd(posts_df[[c_name]] / sd(param_df[[c_name]])))



c_name <- "sfs1_mean"
tibble(
  value = c(param_df[[c_name]], posts_df[[c_name]], 0),
  distribution = c(
    rep("prior", nrow(param_df)), 
    rep("posterior", nrow(posts_df)), 
    "simulated")) %>% 
  ggplot(aes(x = value, colour = distribution)) +
  geom_density(fill=NA, lwd = 2) +
  geom_vline(xintercept = target_df[[c_name]], lwd = 2) +
  scale_colour_manual(values = c("prior" = "grey60", "posterior" = "dodgerblue", "simulated" = "black")) 

den_prior <- density(param_df[[c_name]])
den_post <- density(posts_df[[c_name]])

plot(den_prior, 
     xlim = range(c(range(den_prior$x), range(den_post$x))),
     ylim = range(c(range(den_prior$y), range(den_post$y))),
     col = "red", lwd = 2
)
points(mean(param_df[[c_name]]), 1, pch = "|", cex = 3, col = "red")
lines(den_post, col = "blue", lwd = 2)
points(mean(posts_df[[c_name]]), 1, pch = "|", cex = 3, col = "blue")
points(target_df[[c_name]], 1, pch = "|", cex = 3)
legend()


diff(quantile(param_df[[c_name]], c(0.025, 0.975))) > diff(quantile(posts_df[[c_name]], c(0.025, 0.975)))

quantile(param_df[[c_name]], c(0.025, 0.975))
quantile(posts_df[[c_name]], c(0.025, 0.975))

mean(posts_df[[c_name]] > target_df[[c_name]])
posts_df <- data.frame(res$adj.values)
mean(posts_df[[c_name]] > target_df[[c_name]])

plot(posts_df[["sfs1_mean"]], posts_df[["sfs2_mean"]])
mean(posts_df[["sfs1_mean"]] > posts_df[["sfs2_mean"]])
abline(0,1)

mean(posts_df[["N0"]] > posts_df[["Nb"]])
plot(posts_df[["N0"]], posts_df[["Nb"]])
abline(0,1)

target_df
