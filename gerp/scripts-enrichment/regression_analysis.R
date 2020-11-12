args <- commandArgs(trailingOnly=TRUE)
project_path=args[1]
library(dplyr)
library(data.table)

# read in the genetic map
map <- read.table(file=paste(project_path, "data/annotations/ogut_fifthcM_map_agpv5.bed", sep=""))
colnames(map) <- c("chr", "start", "stop", "s_marker", "m_marker", "useless", "cm")
# read in windows with insertions
windows_ins <- fread(file=paste(project_path, "analyses/windows/windows_ins.bed", sep=""))
windows_ins <- windows_ins %>%
select(V1,V2,V3,V7)
colnames(windows_ins) <- c("chr", "start", "stop", "ins_overlap")
# get the number of insertions that overlap each window
windows_ins <- windows_ins %>%
  mutate(ins_number = ifelse(ins_overlap == 0, 0, 1)) %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))

# reading in windows with GERP
windows_gerp <- fread(file=paste(project_path,"analyses/windows/windows_gerp.bed",sep=""))
windows_gerp <- windows_gerp %>%
select(V1,V2,V3,V12)
colnames(windows_gerp) <- c("chr", "start", "stop", "gerp_overlap")
# get overlap of each window with GERP elements
windows_gerp <- windows_gerp %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
colnames(windows_gerp)[4] <- "gerp_overlap"

windows <- merge(windows_ins,windows_gerp)

# read in windows with open chromatin
windows_open <- fread(file=paste(project_path,"analyses/windows/windows_open.bed",sep=""))
windows_open <- windows_open %>%
select(V1,V2,V3,V7)
colnames(windows_open) <- c("chr", "start", "stop", "open_overlap")
# get overlap of each window with open chromatin
windows_open <- windows_open %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
colnames(windows_open)[4] <- "open_overlap"

windows <- merge(windows,windows_open)

# read in windows with masked base pairs
windows_masked <- fread(file=paste(project_path, "analyses/windows/windows_masked.bed", sep=""))
windows_masked <- windows_masked %>%
select(V1,V2,V3,V7)
colnames(windows_masked) <- c("chr", "start", "stop", "masked_overlap")
# get overlap of each window with masked base pairs
windows_masked <- windows_masked %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
colnames(windows_masked)[4] <- "masked_overlap"

windows <- merge(windows,windows_masked)

windows <- windows %>%
group_by(chr)

map <- map %>%
group_by(chr)

# shift all genetic distances forward by size of most negative distance to make all distances positive
make_positive <- function(x, na.rm = FALSE) (x - min(x))
map <- map %>%
mutate_at(c("cm"), make_positive) %>%
filter(cm == cummax(cm)) %>%
filter(start == cummax(start))

windows_list <- list()

# impute the genetic distance of each window using spline command
for(i in unique(windows$chr)){
  print(i)
  windows_subset <- windows %>%
  filter(chr==i) %>%
  arrange(chr,start)
  print(dim(windows_subset))
  map_subset <- map %>%
  filter(chr==i) %>%
  arrange(chr,start)
  print(dim(map_subset))
  windows_subset$start_cm <-
  spline(x=map_subset$start,y=map_subset$cm,
    xout=windows_subset$start, method = "hyman")$y
  windows_subset$stop_cm <- spline(x=map_subset$start,y=map_subset$cm,
    xout=windows_subset$stop, method = "hyman")$y
  windows_subset <- windows_subset %>%
  mutate(cm=stop_cm-start_cm) %>%
  filter(cm<.15) %>%
  filter(cm>0)
  windows_list[[i]] <- windows_subset
}

windows <- do.call(rbind, windows_list)

# add deletions
windows_del <- fread(file=paste(project_path, "analyses/windows/windows_del.bed", sep=""))
windows_del <- windows_del %>%
select(V1,V2,V3,V7)
colnames(windows_del) <- c("chr", "start", "stop","del_overlap")
windows_del <- windows_del %>%
  mutate(del_number = ifelse(del_overlap == 0, 0, 1)) %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
windows <- merge(x = windows, y = windows_del, all.x = TRUE)

# add inversions
windows_inv <- fread(file=paste(project_path,"analyses/windows/windows_inv.bed", sep = ""))
windows_inv <- windows_inv %>%
select(V1,V2,V3,V7)
colnames(windows_inv) <- c("chr", "start", "stop","inv_overlap")
windows_inv <- windows_inv %>%
  mutate(inv_number = ifelse(inv_overlap == 0, 0, 1)) %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
windows <- merge(x = windows, y = windows_inv, all.x = TRUE)

saveRDS(windows, paste(project_path, "analyses/windows/windows_w_sv_number.RDS", sep=""))

# convert cm (which is in cm/kb) to cm/Mb
data <- data %>%
  mutate(cm_per_mb = cm*100)

# second part of figure: insertions, inversions, deletions, regressed on recombination rate
# make a deta.frame to hold the coefficients from the regression models
betas <- data.frame(SV=NA,Model=NA,Beta=NA, SE=NA)
newdata2 <- data.frame(
  cm_per_mb = rep(seq(from = min(data$cm_per_mb), to =
                           max(data$cm_per_mb), length.out = 100), 3),
  gerp_overlap = (rep(0, each = 100)),
  open_overlap = rep(0, each = 100),
  masked_overlap = rep(0, each = 100))
newdata2 <- newdata2 %>%
  mutate(length=10000)
# regression for insertions
m_ins <- glm(ins_number_sum ~ gerp_overlap + cm_per_mb + open_overlap + masked_overlap +
               open_overlap:gerp_overlap +
               open_overlap:cm_per_mb +
               gerp_overlap:cm_per_mb +
               masked_overlap:gerp_overlap + masked_overlap:cm_per_mb +
               masked_overlap, data = data, family="quasipoisson")
betas[1,] <- c("Insertions", "Full Model", as.numeric(m_ins$coefficients[2]), sqrt(vcov(m_ins))[2,2])

# predict number of insertions in window based on varying recombination rate
newdata2 <- cbind(newdata2, predict(m_ins, newdata2, type = "response", se.fit=TRUE))
newdata2 <- within(newdata2, {
  Insertion <- (fit)
  Insertion.LL <- (fit - 1.96 * se.fit)
  Insertion.UL <- (fit + 1.96 * se.fit)
})
newdata2 <- newdata2 %>%
  dplyr::select(-c(fit,se.fit,residual.scale))
# regression for inversions
m_inv <- glm(inv_number_sum ~ gerp_overlap + cm_per_mb + open_overlap + masked_overlap +
               open_overlap:gerp_overlap +
               open_overlap:cm_per_mb +
               gerp_overlap:cm_per_mb +
               masked_overlap:gerp_overlap + masked_overlap:cm_per_mb +
               masked_overlap, data = data, family="quasipoisson")
betas[3,] <- c("Inversions", "Full Model", as.numeric(m_inv$coefficients[2]), sqrt(vcov(m_inv))[2,2])
# prediction for inversions
newdata2 <- cbind(newdata2, predict(m_inv, newdata2, type = "response", se.fit=TRUE))
newdata2 <- within(newdata2, {
  Inversion <- (fit)
  Inversion.LL <- (fit - 1.96 * se.fit)
  Inversion.UL <- (fit + 1.96 * se.fit)
})
newdata2 <- newdata2 %>%
  dplyr::select(-c(fit,se.fit,residual.scale))
# regression for deletions
m_del <- glm(del_number_sum ~ gerp_overlap + cm_per_mb + open_overlap + masked_overlap +
               open_overlap:gerp_overlap +
               open_overlap:cm_per_mb +
               gerp_overlap:cm_per_mb +
               masked_overlap:gerp_overlap + masked_overlap:cm_per_mb +
               masked_overlap, data = data, family="quasipoisson")
betas[2,] <- c("Deletions", "Full Model", as.numeric(m_del$coefficients[2]), sqrt(vcov(m_del))[2,2])

# prediction for deletions
newdata2 <- cbind(newdata2, predict(m_del, newdata2, type = "response", se.fit=TRUE))
newdata2 <- within(newdata2, {
  Deletion <- (fit)
  Deletion.LL <- (fit - 1.96 * se.fit)
  Deletion.UL <- (fit + 1.96 * se.fit)
})

# making the figures

library(RColorBrewer)
library(tidyr)
library(ggplot2)

# multiplot code from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# plotting the coefficients of the regressions models
betas$Beta <- as.numeric(betas$Beta)
betas$SE <- as.numeric(betas$SE)
betas_filtered <- betas %>%
  filter(Model == "Full Model")
p1 <- ggplot(betas_filtered, aes(x=SV, y=Beta, colour=Model)) +
    geom_errorbar(aes(ymin=Beta-SE, ymax=Beta+SE), width=.1) +
    geom_point() +
  labs(tag="A") +
    geom_hline(yintercept = 0, linetype = "dashed")   +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        legend.position = "none")

# Plotting the prediction for each type of SV by varying recombination rate
newdata2_ins <- newdata2 %>%
  dplyr::select(cm_per_mb,Insertion,Insertion.UL,Insertion.LL)
newdata2_ins$SV <- rep("Insertion", length(newdata2_ins$cm_per_mb))
newdata2_ins <- newdata2_ins %>%
  rename(Estimate = Insertion, UL = Insertion.UL, LL = Insertion.LL)
newdata2_del <- newdata2 %>%
  dplyr::select(cm_per_mb,Deletion,Deletion.UL,Deletion.LL)
newdata2_del$SV <- rep("Deletion", length(newdata2_del$cm_per_mb))
newdata2_del <- newdata2_del %>%
  rename(Estimate = Deletion, UL = Deletion.UL, LL = Deletion.LL)
newdata2_inv <- newdata2 %>%
  dplyr::select(cm_per_mb,Inversion,Inversion.UL,Inversion.LL)
newdata2_inv$SV <- rep("Inversion", length(newdata2_inv$cm_per_mb))
newdata2_inv <- newdata2_inv %>%
  rename(Estimate = Inversion, UL = Inversion.UL, LL = Inversion.LL)
newdata2_long <- rbind(newdata2_del,newdata2_inv,newdata2_ins)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p2 <- ggplot(newdata2_long, aes(cm_per_mb, Estimate, fill = SV, color =SV)) +
  geom_line(size = .5) +
  labs(x = "cM per Mb", y = "Number overlapping", tag="B") +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25)  +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold")) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme(legend.title=element_text(size=14),
    legend.text=element_text(size=13))

multiplot(p1,p2)
