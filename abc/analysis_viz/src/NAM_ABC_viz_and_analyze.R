library(tidyverse)
library(patchwork)
library(cowplot)
library(broom)
library(lme4)
theme_set(theme_cowplot(font_size = 22))
select <- dplyr::select

params <- 
  c("Na", "N0", "Nb", "B_t", "mu_sv",
    "sfs1_shape", "sfs1_mean", "sfs2_shape", "sfs2_mean")


param_df <- read_delim(
  "predict/pop.txt",
  col_names = params,
  delim = " "
) %>% 
  dplyr::select(-c(Na, B_t))


range(param_df$N0) / 1000
unique(param_df$mu_sv)


all_abc_post <-
    read_delim(file = "predict/data/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_df.csv", delim = ",") %>% 
  mutate(s_diff = sfs2_mean - sfs1_mean)


all_abc_post %>% 
  group_by(mu_sv) %>% 
  summarise(m1 = mean(sfs1_mean),
            m2 = mean(sfs2_mean), 
            mean(s_diff)) %>% 
  mutate(m2 - m1)



1 - mean(all_abc_post$sfs1_mean > 0.05)
1 - mean(all_abc_post$sfs2_mean > 0.05)



alph <- 0.95
low <- (1 - alph) / 2
high <- 1 - low

cred_df <- 
  all_abc_post %>% 
  select(sfs1_mean, sfs2_mean, window, mu_sv) %>%
  rename(SNPs = sfs1_mean, SVs = sfs2_mean) %>% 
  pivot_longer(cols = c(SNPs, SVs), 
               names_to = "mutation_type", 
               values_to = "s") %>% 
  mutate(s = -s) %>% 
  #group_by(window, mutation_type, mu_sv) %>% 
  group_by(window, mutation_type) %>% 
  summarise(qlow = quantile(s, low),
            qhigh = quantile(s, high),
            qmed = quantile(s, 0.5)) %>%
  mutate(x = rnorm(n(), 0, 0.1) + as.numeric(as.factor(mutation_type)),
         mx = as.numeric(as.factor(mutation_type))) %>%
  group_by(mutation_type) %>% 
  mutate(mean = mean(qmed))


cred_df %>% 
  ggplot(aes(mutation_type, qmed)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.2) +
  #facet_wrap(~mu_sv) +
  ylab("s") +
  xlab("mutation type") +
  ggsave("analysis_viz/figures/post_s_fixedMu.pdf", width = 5, height = 6)
  



win_df <- read_delim(
  "predict/data/NAM_abc_out/window_stats_20Mb.txt", delim = "\t", skip = 1,
  col_names = c("chrom", "start", "end", "fold0", "fold4", "cds_bps", "genes")) %>% 
  mutate(window = 1:n(),
         length = end - start)



summarise_if(param_df, .predicate = is.numeric, .funs = mean) %>% 
  pivot_longer(cols = everything(), names_to = "param", values_to = "mean_prior")




sfs1 <- 
  all_abc_post %>%
  ggplot(aes(x = -sfs1_mean, group = window)) +
  geom_density(colour=alpha("black", 0.2)) +
  geom_density(aes(x = -sfs1_mean), inherit.aes = F) + 
  geom_density(data = param_df, aes(x = -sfs1_mean), lwd = 1.1, inherit.aes = F, colour = "red") +
  xlim(-0.1, 0.01) +
  ylim(0, 300) +
  xlab("selection 0-fold")
sfs1

sfs2 <- 
  all_abc_post %>%
  ggplot(aes(x = -sfs2_mean, group = window)) +
  geom_density(colour=alpha("black", 0.2)) +
  #geom_density(aes(x = -sfs1_mean), inherit.aes = F) +
  geom_density(aes(x = -sfs2_mean), colour = "blue", inherit.aes = F) + 
  geom_density(data = param_df, aes(x = -sfs2_mean), lwd = 1.1, inherit.aes = F, colour = "red") +
  xlim(-0.1, 0.1) +
  ylim(0, 300) +
  xlab("selection SV")
sfs2

all_abc_post %>%
  filter(window == sample(1:113, 1)) %>% 
  ggplot(aes(x = -sfs2_mean, group = window)) +
  geom_density(colour= alpha("black", 0.2)) +
  #geom_density(aes(x = sfs1_mean), inherit.aes = F) +
  #geom_density(aes(x = sfs2_mean), colour = "blue", inherit.aes = F) + 
  geom_density(data = param_df, aes(x = -sfs2_mean), lwd = 1.1, inherit.aes = F, colour = "red") +
  xlim(-0.1, 0.01) +
  xlab("selection SV")



all_abc_post %>%
  filter(window == sample(1:113, 1)) %>% 
  ggplot(aes(x = -sfs1_mean, group = window)) +
  geom_density(colour= alpha("black", 0.2)) +
  #geom_density(aes(x = sfs1_mean), inherit.aes = F) +
  #geom_density(aes(x = sfs2_mean), colour = "blue", inherit.aes = F) + 
  geom_density(data = param_df, aes(x = -sfs1_mean), lwd = 1.1, inherit.aes = F, colour = "red") +
  xlim(-0.1, 0.01) 




N0 <- 
  all_abc_post %>%
  ggplot(aes(x = log10(N0), group = window)) +
  geom_density(colour=alpha("black", 0.2)) +
  #geom_density(aes(x = N0), inherit.aes = F) + 
  geom_density(data = param_df, aes(x = log10(N0)), lwd = 1.1, inherit.aes = F, colour = "red") +
  xlim(2, 5)
N0

Nb <- 
  all_abc_post %>%
  ggplot(aes(x = Nb, group = window)) +
  geom_density(colour=alpha("black", 0.2)) +
  #geom_density(aes(x = Nb), inherit.aes = F) + 
  geom_density(data = param_df, aes(x = Nb), lwd = 1.1, inherit.aes = F, colour = "red") +
  ylim(0, 0.01) +
  xlim(0, 2000)
Nb

mu <- 
  all_abc_post %>%
  ggplot(aes(x = mu_sv, group = window)) +
  geom_density(colour=alpha("black", 0.2)) +
  geom_density(aes(x = mu_sv), inherit.aes = F) + 
  geom_density(data = param_df, aes(x = mu_sv), lwd = 1.1,  inherit.aes = F, colour = "red") +
  xlab(expression(mu[sv]))
  #xlim(0, 5e-8)
mu


sv_sums  = c(252,239,188,166,149,252,359,330,230,151,210,287,298,256,181,125,52,115,149,92,106,125,64,58,191,189,127,182,254,246,353,175,196,122,146,200,40,137,123,213,285,296,226,157,139,129,244,228,183,137,134,213,193,300,226,208,122,134,175,245,203,110,132,174,185,156,162,257,149,149,114,130,126,65,46,81,161,139,112,84,202,181,143,19,80,94,211,261,141,98,97,105,120,55,70,76,219,149,114,109,117,120,100,17,96,91,106,110,120,102,108,108,23)
sv_counts <- 
  tibble(
    sv_counts = sv_sums,
    window = seq_along(sv_sums)
  )



post_window_mean <-
  all_abc_post %>%
  group_by(window, sfs) %>%
  mutate(pr_dfe2gtdfe1 = sfs1_mean < sfs2_mean) %>%
  summarise_all(.funs = mean) %>% 
  mutate(idx = 1:n()) %>% 
  full_join(., win_df, by = "window") %>% 
  filter(length >= 20000000) %>% 
  left_join(. , sv_counts, by = "window")




post_window_mean %>% 
  ggplot(aes(mu_sv)) +
  geom_histogram()


wilcox.test(post_window_mean$sfs1_mean, post_window_mean$sfs2_mean, paired = T)
mean(post_window_mean$sfs1_mean)
mean(post_window_mean$sfs2_mean)
mean(post_window_mean$sfs1_mean - post_window_mean$sfs2_mean)


full_join(
  post_window_mean %>%
    ungroup() %>% 
    select(N0, Nb, mu_sv, sfs1_mean, sfs2_mean, sfs1_shape, sfs2_shape) %>% 
    summarise_all(mean) %>% 
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "mean"), 
  post_window_mean %>%
    ungroup() %>% 
    select(N0, Nb, mu_sv, sfs1_mean, sfs2_mean, sfs1_shape, sfs2_shape) %>% 
    summarise_all(sd) %>% 
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "sd")
) %>% 
  xtable::xtable(display = c("f", "f", "f", "f"), digits = 4) %>% 
  xtable::print.xtable(include.rownames = FALSE)


a <- post_window_mean %>%
  mutate(dfe_0fold = sfs1_mean, 
         dfe_sv = sfs2_mean) %>%
  #filter(mu_sv > 1e-8) %>% 
  dplyr::select(sfs, dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "s") %>% 
  ggplot(aes(dfe, -s, fill = dfe)) +
  ylab("mean posterior s") +
  xlab("Variant type") +
  scale_x_discrete(labels = c('0-fold','SV')) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
  theme(legend.position = "n") +
  scale_fill_manual(values = c("dfe_0fold" = "dodgerblue", "dfe_sv" = "forestgreen"))
a



b <- post_window_mean %>%
  mutate(dfe_0fold = sfs1_mean, 
         dfe_sv = sfs2_mean) %>% 
  dplyr::select(sfs, dfe_0fold, dfe_sv) %>% 
  mutate(diff = dfe_sv -dfe_0fold) %>% 
  ggplot() +
  geom_histogram(
    aes(diff), bins = 10,
    fill = "grey90", colour = "black"
  ) +
  xlim(-0.1, 0.1) +
  xlab(expression(paste(s[sv] - s["0fold"])))
b



post_lm_df <- 
  post_window_mean %>% 
  pivot_longer(cols = c(sfs2_mean, sfs1_mean), 
               names_to = "mutation_type", values_to = "s") %>%
  mutate(cds_bps = cds_bps/1e5)

post_lm_df %>% 
  mutate(
    mutation_type = str_replace(mutation_type, "sfs1_mean", "s_0"),
    mutation_type = str_replace(mutation_type, "sfs2_mean", "s_sv")
    ) %>% 
  ggplot(aes(cds_bps, -s, colour = mutation_type)) +
  geom_point() +
  geom_smooth(se = F, method = "lm") +
  scale_colour_manual(
    values = c("s_0" = "dodgerblue", "s_sv" = "forestgreen")
  ) +
  xlab(expression(paste("CDS base pairs (×", 10^5, ")"))) +
  ylab("s") +
  ggsave(filename = "analysis_viz/figures/mut_cds.pdf", width = 10, height = 10)


m_full <- lm(-s ~ cds_bps*mutation_type + cds_bps + mutation_type, data = post_lm_df)
summary(m_full)
m_fix <- lm(-s ~ cds_bps + mutation_type, data = post_lm_df)
summary(m_fix)
m_simp <- lm(-s ~ cds_bps, data = post_lm_df)
summary(m_simp)
m_int <- lm(-s ~ 1, data = post_lm_df)
summary(m_int)

anova(m_simp, m_int, test="LRT")
anova(m_fix, m_simp, test="LRT")
anova(m_full, m_fix, test="LRT")
