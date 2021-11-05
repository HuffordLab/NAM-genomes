library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
summary.svs.filename <- args[1]
summary.rils.filename <- args[2]
sv.filter <- as.numeric(args[3])
ril.filter <- as.numeric(args[4])

# setwd("~/projects/sv_nams/analysis/projection/sv_filters")
# # setwd("analysis/projection/sv_filters")
# summary.svs.filename <- "tassel_summary_NAM_rils_projected_svs_10-perc-missing-filter3.txt"
# summary.rils.filename <- "tassel_summary_NAM_rils_projected_svs_10-perc-missing-filter4.txt"
# sv.filter <- 10
# ril.filter <- 5

# plot distribution svs
summary.svs <- fread(summary.svs.filename, header = TRUE, data.table = FALSE)

plot.missing.svs <- ggplot(summary.svs, aes(x = `Proportion Missing`)) +
  geom_histogram(binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Proportion missing",
       y = "Count (SVs)",
       title = paste0(NROW(summary.svs), " SVs with < ", sv.filter, "% missing data"))

outfile.sv.dist <- paste0("SVs-", sv.filter, "-perc-missing_RILs-", ril.filter, "-perc-missing_distribution-SVs.png")
ggsave(plot = plot.missing.svs, filename = outfile.sv.dist, device = "png")


# plot distribution rils
summary.rils <- fread(summary.rils.filename, header = TRUE, data.table = FALSE)

plot.missing.rils <- ggplot(summary.rils, aes(x = `Proportion Missing`)) +
  geom_histogram(binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Proportion missing",
       y = "Count (RILs)",
       title = paste0(sum(summary.rils[, "Proportion Missing"] < ril.filter/100), " RILs with < ", ril.filter, "% missing data"),
       subtitle = paste0("(SVs with < ", sv.filter, "% missing data)"))

outfile.ril.dist <- paste0("SVs-", sv.filter, "-perc-missing_RILs-", ril.filter, "-perc-missing_distribution-RILs.png")
ggsave(plot = plot.missing.rils, filename = outfile.ril.dist, device = "png")


# plot rils kept per family
rils.kept <- summary.rils[which(summary.rils[, "Proportion Missing"] < ril.filter/100), "Taxa Name"]
rils.per.family <- data.frame(table(substr(rils.kept, start = 1, stop = 4)))
rils.per.family$Var1 <- factor(rils.per.family$Var1,
                               levels = c("Z001", "Z002", "Z003", "Z004", "Z005", "Z006", "Z007",
                                          "Z008", "Z009", "Z010", "Z011", "Z012", "Z013", "Z014",
                                          "Z015", "Z016", "Z018", "Z019", "Z020", "Z021", "Z022",
                                          "Z023", "Z024", "Z025", "Z026"))

plot.rils.per.family <- ggplot(rils.per.family, aes(x = Var1, y = Freq)) + 
  geom_col() +
  scale_x_discrete(drop = FALSE) +
  labs(x = "NAM families",
       y = "RILs",
       title = "RILs kept per family after filtering",
       subtitle = paste0("(SVs < ", sv.filter, "% missing, RILs < ", ril.filter, "% missing)"))

outfile.ril.kept <- paste0("SVs-", sv.filter, "-perc-missing_RILs-", ril.filter, "-perc-missing_RILs-kept.png")
ggsave(plot = plot.rils.per.family, filename = outfile.ril.kept, device = "png", width = 9.5)

