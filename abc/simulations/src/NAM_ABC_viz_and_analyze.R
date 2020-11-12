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
  #"NAM_DFE_ABC/pop_fixedMu.txt",
  "NAM_DFE_ABC/pop.txt",
  col_names = params,
  delim = " "
) %>% 
  dplyr::select(-c(Na, B_t))


range(param_df$N0) / 1000
unique(param_df$mu_sv)
# 
# all_abc_post <-
#   bind_rows(
#     read_delim(file = "NAM_DFE_ABC/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_fixedMu_1e-08_df.csv", delim = ",") %>% 
#       mutate(mu_sv = "1e-08"),
#     read_delim(file = "NAM_DFE_ABC/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_fixedMu_1e-09_df.csv", delim = ",")%>% 
#       mutate(mu_sv = "1e-09"),
#     read_delim(file = "NAM_DFE_ABC/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_fixedMu_1e-10_df.csv", delim = ",")%>% 
#       mutate(mu_sv = "1e-10")
#   ) %>% 
#   mutate(s_diff = sfs2_mean - sfs1_mean)


all_abc_post <-
    read_delim(file = "NAM_DFE_ABC/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_df.csv", delim = ",") %>% 
  mutate(s_diff = sfs2_mean - sfs1_mean)

# full_join(
# all_abc_post %>%
#   select(-window, -s_diff) %>% 
#   summarise_if(is.numeric, .funs = mean) %>% 
#   pivot_longer(cols = everything(), names_to = "parameter", values_to = "mean"), 
# all_abc_post %>%
#   select(-window, -s_diff) %>%
#   summarise_if(is.numeric, .funs = sd) %>% 
#   pivot_longer(cols = everything(), names_to = "parameter", values_to = "sd") 
# ) %>% 
#   xtable::xtable(display = c("f", "f", "f", "f"), digits = 4) %>% 
#   xtable::print.xtable(include.rownames = FALSE)

all_abc_post %>% 
  group_by(mu_sv) %>% 
  summarise(m1 = mean(sfs1_mean),
            m2 = mean(sfs2_mean), 
            mean(s_diff)) %>% 
  mutate(m2 - m1)



1 - mean(all_abc_post$sfs1_mean > 0.05)
1 - mean(all_abc_post$sfs2_mean > 0.05)


# c("1e-08", "1e-09", "1e-10") %>% 
#   map(~{
#     ndf <- 
#       all_abc_post %>%
#       group_by(window, mu_sv) %>% 
#       summarise(sfs1_mean = mean(sfs1_mean), sfs2_mean = mean(sfs2_mean)) %>% 
#       filter(mu_sv == .x) 
#     c(mean(ndf$sfs1_mean) - mean(ndf$sfs2_mean), wilcox.test(ndf$sfs1_mean, ndf$sfs2_mean, paired = T)$p.value)
#     #wilcox.test(ndf$sfs1_mean, ndf$sfs2_mean, paired = T)
#   })
  

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
  ggsave("NAM_DFE_ABC/figures/post_s_fixedMu.pdf", width = 5, height = 6)
  

# #wrong
# cred_df %>% 
#   ggplot() +
#   geom_segment(aes(x = x, xend = x, y = qlow, yend = qhigh), alpha = 0.2) +
#   geom_point(aes(x, qmed, colour = mutation_type)) +
#   geom_point(data = distinct(select(cred_df, mx, mean)), 
#             aes(mx, mean, fill = mutation_type), size = 4, shape = 21) +
#   scale_fill_manual(values = c("SNPs" = "dodgerblue", "SVs" = "forestgreen")) +
#   scale_colour_manual(values = c("SNPs" = "dodgerblue", "SVs" = "forestgreen")) +
#   ylab("s") +
#   facet_grid(~mu_sv) +
#   scale_x_continuous(name = "",
#                      breaks = c(1, 2),
#                      labels=c("SNPs", "SVs"))+
#   theme(legend.position = "n") +
#   ggsave("NAM_DFE_ABC/figures/post_s.png")
  
  
cred_df %>% 
  mutate(interval = 1/(qhigh - qlow)) %>% 
  ggplot(aes(mutation_type, qmed, size = interval, colour = mutation_type)) + 
  ggforce::geom_sina() +
  scale_fill_manual(values = c("SNPs" = "dodgerblue", "SVs" = "forestgreen")) +
  scale_colour_manual(values = c("SNPs" = "dodgerblue", "SVs" = "forestgreen")) +
  ylab("s") +
  theme(legend.position = "n") +
  ggsave("NAM_DFE_ABC/figures/post_s.png")
  


win_df <- read_delim(
  "NAM_DFE_ABC/NAM_abc_out/window_stats_20Mb.txt", delim = "\t", skip = 1,
  col_names = c("chrom", "start", "end", "fold0", "fold4", "cds_bps", "genes")) %>% 
  mutate(window = 1:n(),
         length = end - start)


# sfs <- read_delim("~/Dropbox/NAM_DFE_ABC/abc_out/NAM_SFS_INS_DEL-gte2Kb_all.txt", "\t", col_names = FALSE)
# sums <- apply(sfs, 1, sum)



pgtep <- 
  full_join(
    pivot_longer(all_abc_post, cols = -c(window, sfs), names_to = "param", values_to = "post"),
    summarise_if(param_df, .predicate = is.numeric, .funs = mean) %>% 
      pivot_longer(cols = everything(), names_to = "param", values_to = "mean_prior"),
    by = "param") %>%
  group_by(window, param) %>% 
  summarise(post_gte_pre = mean(post >= mean_prior)) %>% 
  ggplot(aes(post_gte_pre)) +
  geom_histogram() +
  facet_wrap(~param, scales = "free", nrow = 2) +
  xlab(expression(paste(posterior[i], " ≥ ", bar(prior))))


param_sd <- 
  full_join(
  all_abc_post %>% 
    group_by(window) %>% 
    summarise_if(.predicate = is.numeric,.funs = sd) %>% 
    pivot_longer(cols = -window, names_to = "param", values_to = "sd_post"),
  summarise_if(param_df, .predicate = is.numeric,.funs = sd) %>% 
    pivot_longer(cols = everything(), names_to = "param", values_to = "sd_prior"),
  by = "param")

param_mean <- 
  full_join(
  all_abc_post %>% 
    group_by(window) %>% 
    summarise_if(.predicate = is.numeric,.funs = mean) %>% 
    pivot_longer(cols = -window, names_to = "param", values_to = "mean_post"),
  summarise_if(param_df, .predicate = is.numeric, .funs = mean) %>% 
    pivot_longer(cols = everything(), names_to = "param", values_to = "mean_prior"),
  by = "param")


full_join(
  param_sd %>% 
    group_by(param, window) %>% 
    summarise(lt = mean(sd_post < sd_prior),
              diff = mean(sd_prior - sd_post) /sd_prior) %>% 
    ungroup() %>% 
    group_by(param) %>% 
    summarise(sd_prop_less = mean(lt == 1), 
              sd_mean_diff = mean(diff)),
  
  param_mean %>% 
    group_by(param, window) %>% 
    summarise(lt = mean(mean_post < mean_prior),
              diff = mean(mean_prior - mean_post) / mean_prior) %>% 
    ungroup() %>% 
    group_by(param) %>% 
    summarise(mean_prop_less = mean(lt == 1), 
              mean_mean_diff = mean(diff)),
  by = "param"
  
)


summarise_if(param_df, .predicate = is.numeric, .funs = mean) %>% 
  pivot_longer(cols = everything(), names_to = "param", values_to = "mean_prior")


xx <- seq(0, quantile(all_abc_post$sfs2_mean, 0.99), length.out = 500)
DFE2_vec <- sample(1:nrow(all_abc_post), 1) %>% 
  map(~{
    dgamma(x = xx, shape = all_abc_post$sfs2_shape[.x], scale = all_abc_post$sfs2_mean[.x]/all_abc_post$sfs2_shape[.x])
  })
mxx <- max(map_dbl(DFE2_vec, max))
plot(xx, xx, type = 'n', ylim=c(0, 1))
walk(DFE2_vec, ~{lines(xx, .x/max(.x), col = alpha("black", 0.9))})

xx <- seq(0, quantile(all_abc_post$sfs1_mean, 0.99), length.out = 500)
DFE2_vec <- sample(1:nrow(all_abc_post), 1) %>% 
  map(~{
    dgamma(x = xx, shape = all_abc_post$sfs1_shape[.x], scale = all_abc_post$sfs1_mean[.x]/all_abc_post$sfs1_shape[.x])
  })
mxx <- max(map_dbl(DFE2_vec, max))
#plot(xx, xx, type = 'n', ylim=c(0, mxx))
walk(DFE2_vec, ~{lines(xx, .x/max(.x), col = alpha("blue", 0.9))})


wilcox.test(all_abc_post$sfs1_mean, all_abc_post$sfs2_mean, paired = TRUE)
mean(all_abc_post$sfs1_mean)
mean(all_abc_post$sfs2_mean)
mean(all_abc_post$sfs1_mean < all_abc_post$sfs2_mean)


all_abc_post %>% 
  dplyr::select(sfs1_mean, sfs2_mean) %>% 
  pivot_longer(cols = everything(), values_to = "s", names_to = "sfs") %>% 
  ggplot() +
  geom_boxplot(aes(sfs, -s))

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

hist(all_abc_post$mu_sv)

all_abc_post %>%
  filter(window %in% sample(1:103, 15, replace = F)) %>% 
  ggplot(aes(mu_sv, sfs2_mean, colour = factor(window))) +
  geom_point( alpha = 1) +
  geom_smooth(method = 'lm', se = F, colour = "black") +
  facet_wrap(~window) +
  theme(legend.position = 'n') +
  ggsave("NAM_DFE_ABC/figures/post_musv_SVs.png")


slope_df <- all_abc_post %>%
  group_by(window) %>% 
  do({
    mod <- lm(-sfs2_mean ~ mu_sv, data = .)
    slope <- coef(mod)[2]
    tibble(slope = slope)
  })

slope_df %>% 
  ggplot(aes(scale(slope))) +
  geom_histogram() +
  ggsave("NAM_DFE_ABC/figures/muSV_slope_SVs.png")



sfs1 + sfs2 + plot_layout(nrow = 2) + ggsave("NAM_DFE_ABC/figures/post_sel.png")

mean(all_abc_post$sfs2_mean > all_abc_post$sfs1_mean)


all_abc_post %>% 
  ggplot(aes(sfs2_mean, sfs1_mean, colour = window)) +
  geom_point(alpha = 0.07) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_smooth() +
  theme(legend.position = "n")


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
  ggplot(aes(mu_sv, sfs2_mean - sfs1_mean)) +
  geom_point() +
  geom_smooth()

post_window_mean %>% 
  ggplot(aes(mu_sv, sfs2_mean)) +
  geom_point() +
  geom_smooth()


post_window_mean %>%
  ggplot(aes(mu_sv, sfs1_mean)) +
  geom_point() +
  geom_smooth()


mean(post_window_mean$mu_sv) / 100

summary(lm(post_window_mean$sfs2_mean ~ post_window_mean$mu_sv))
summary(lm(post_window_mean$sfs1_mean ~ post_window_mean$mu_sv))

post_window_mean %>% 
  ggplot(aes(mu_sv)) +
  geom_histogram()


post_window_mean %>%
  mutate(dfe_0fold = sfs1_mean*1000, 
         dfe_sv = sfs2_mean*1000) %>% 
  dplyr::select(sfs, dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "Ns") %>% 
  ggplot(aes(Ns, dfe)) +
  ggridges::geom_density_ridges(alpha = 0.5, rel_min_height = 0.001) +
  ggsave(filename = "~/Desktop/dfe.png")



t.test(filter(post_window_mean, mu_sv > 1e-8)$sfs1_mean, filter(post_window_mean, mu_sv > 1e-8)$sfs2_mean, paired = T, var.equal = T)
wilcox.test(filter(post_window_mean, mu_sv > 1e-8)$sfs1_mean, filter(post_window_mean, mu_sv > 1e-8)$sfs2_mean)
mean(filter(post_window_mean, mu_sv > 1e-8)$sfs1_mean)
mean(filter(post_window_mean, mu_sv > 1e-8)$sfs2_mean)
mean(filter(post_window_mean, mu_sv > 1e-8)$sfs1_mean - filter(post_window_mean, mu_sv > 1e-8)$sfs2_mean)


t.test(post_window_mean$sfs1_mean, post_window_mean$sfs2_mean, paired = T, var.equal = T)
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

post_window_mean %>%
  mutate(dfe_0fold = sfs1_mean, 
         dfe_sv = sfs2_mean) %>%
  dplyr::select(sfs, dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "s") %>% 
  ggplot(aes(dfe, -s, colour = dfe)) +
  ggforce::geom_sina() +
  scale_colour_manual(
    values = c("dfe_0fold" = "dodgerblue", "dfe_sv" = "forestgreen")
    )

plot(post_window_mean$sfs1_shape, post_window_mean$sfs2_shape)

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

# a + b +  plot_layout(nrow = 2) + ggsave("NAM_DFE_ABC/figures/dfe.pdf", width = 7, height = 12)
# a + b +  plot_layout(nrow = 1) + ggsave("NAM_DFE_ABC/figures/dfe.png", width = 12, height = 6)


all_abc_post %>%
  mutate(dfe_0fold = sfs1_mean, 
         dfe_sv = sfs2_mean) %>% 
  dplyr::select(dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "Ns") %>% 
  ggplot(aes(Ns, dfe)) +
  ggridges::geom_density_ridges(alpha = 0.5, rel_min_height = 0.001)

post_window_mean %>%
  mutate(dfe_0fold = sfs1_mean, 
         dfe_sv = sfs2_mean) %>% 
  dplyr::select(sfs, dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "s") %>% 
  ggplot(aes(-s, dfe)) +
  ggridges::geom_density_ridges(alpha = 0.5, rel_min_height = 0.001) +
  ggsave(filename = "~/Desktop/dfe.png")

post_window_mean %>%
  ggplot(aes(sfs1_mean, sfs2_mean, colour = sfs)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

post_window_mean %>%
  ggplot(aes(window, sfs2_mean, colour = chrom)) +
  geom_point() +
  ylab("s") +
  ggsave(filename = "~/Desktop/s_window.png")

post_window_mean %>%
  ggplot(aes(window, sfs1_mean, colour = chrom)) +
  geom_point() +
  ggsave(filename = "~/Desktop/s_window.png")

post_window_mean %>%
  ggplot(aes(genes, -sfs1_mean)) +
  geom_point(aes(colour = cds_bps)) +
  geom_smooth(se = F, colour = "grey70") +
  ylab("0-fold s") +
  xlab("number of genes") +
  ggsave(filename = "~/Desktop/fold0_gene_window.png")

post_window_mean %>%
  ggplot(aes(cds_bps, -sfs2_mean)) +
  geom_point() +
  geom_smooth(se = F, colour = "grey70") +
  ylab("SV s") +
  xlab("cds bps") +
  ggsave(filename = "NAM_DFE_ABC/figures/SV_gene_window.png")


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
  #facet_wrap(~mu_sv) +
  #xlab(expression(paste0("CDS base pairs (", 1*10^5)) +
  xlab(expression(paste("CDS base pairs (×", 10^5, ")"))) +
  ylab("s") +
  ggsave(filename = "NAM_DFE_ABC/figures/mut_cds.pdf", width = 10, height = 10)

# post_window_mean %>% 
#   select(cds_bps, sfs2_mean, sfs1_mean, mu_sv) %>% 
#   pivot_longer(cols = c(sfs2_mean, sfs1_mean), names_to = "mutation_type", values_to = "s") %>% 
#   #group_by(mutation_type, mu_sv) %>%
#   group_by(mutation_type) %>%
#   group_modify(~{
#     lm(-s ~ cds_bps, data = .) %>% broom::tidy()
#   }) %>% 
#   filter(term != "(Intercept)") %>% 
#   select(-std.error, -statistic) %>% 
#   xtable::xtable()


m_full <- lm(-s ~ cds_bps*mutation_type + cds_bps + mutation_type, data = post_lm_df)
summary(m_full)
m_fix <- lm(-s ~ cds_bps + mutation_type, data = post_lm_df)
summary(m_fix)
m_simp <- lm(-s ~ cds_bps, data = post_lm_df)
summary(m_simp)
m_int <- lm(-s ~ 1, data = post_lm_df)
summary(m_int)


plot(post_lm_df$cds_bps, predict(m_full))
plot(post_lm_df$cds_bps, predict(m_fix))
plot(post_lm_df$cds_bps, predict(m_simp))
plot(post_lm_df$cds_bps, predict(m_int))

anova(m_simp, m_int, test="LRT")
anova(m_fix, m_simp, test="LRT")
anova(m_full, m_fix, test="LRT")

post_window_mean %>%
  ggplot(aes(sv_counts, -sfs2_mean, colour = cds_bps)) +
  geom_point() +
  geom_smooth(se = F, colour = "grey70") +
  ggsave("NAM_DFE_ABC/figures/svDensity_vs_s.png")


summary(lm(-sfs2_mean ~ sv_counts + cds_bps, data = post_window_mean))

post_window_mean %>%
  ggplot(aes(cds_bps, -sfs1_mean, colour = cds_bps)) +
  geom_point() +
  geom_smooth(se = F, colour = "grey70")

post_window_mean %>%
  ggplot(aes(cds_bps, sv_counts)) +
  geom_point() +
  geom_smooth(se = T, method = "lm", colour = "grey70") +
  ggsave("NAM_DFE_ABC/figures/svCount_vs_cds.png")


summary(lm(sfs2_mean ~ cds_bps, data = post_window_mean))
summary(lm(sfs1_mean ~ cds_bps, data = post_window_mean))


mod_full <- post_window_mean %>% 
  select(cds_bps, sfs2_mean, sfs1_mean) %>% 
  pivot_longer(cols = c(sfs2_mean, sfs1_mean), names_to = "mutation_type", values_to = "s") %>% 
  lm(-s ~ cds_bps*mutation_type, data = .) %>% 
  summary()

post_window_mean %>% 
  select(cds_bps, sfs2_mean, sfs1_mean) %>% 
  pivot_longer(cols = c(sfs2_mean, sfs1_mean), names_to = "mutation_type", values_to = "s") %>% 
  lm(-s ~ cds_bps + mutation_type, data = .) %>% 
  summary()


length(unique(post_window_mean$window))
mean(post_window_mean$pr_dfe2gtdfe1 >= 0.95)
hist(post_window_mean$pr_dfe2gtdfe1, breaks = 100)
1/mean(1/post_window_mean$pr_dfe2gtdfe1)
exp(mean(log(post_window_mean$pr_dfe2gtdfe1)))
mean(post_window_mean$pr_dfe2gtdfe1)
hist(log10(post_window_mean$sfs1_mean))
hist(log10(post_window_mean$sfs2_mean))
hist(post_window_mean$sfs2_mean)
mean(post_window_mean$sfs2_mean)
mean(post_window_mean$sfs1_mean)
plot(post_window_mean$sfs1_mean, post_window_mean$sfs2_mean)
abline(0,1)
mean(post_window_mean$sfs1_mean < post_window_mean$sfs2_mean)
plot(post_window_mean$Nb, post_window_mean$N0)
hist(log10(post_window_mean$N0), breaks = 100)
hist(log10(post_window_mean$Nb), breaks = 100)
plot(log10(post_window_mean$N0), log10(post_window_mean$Nb))

hist(log10(post_window_mean$N0), breaks = 100)
hist(log10(post_window_mean$Nb), breaks = 100)
plot(post_window_mean$N0, post_window_mean$Nb)
abline(0,1)

plot(post_window_mean$sfs2_mean)
plot(post_window_mean$sfs1_mean)

plot(post_window_mean$sfs2_mean, post_window_mean$sfs1_mean)
abline(0,1)

plot(post_window_mean$sfs2_mean, post_window_mean$mu_sv)

