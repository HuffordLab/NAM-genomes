args <- commandArgs(trailingOnly = TRUE)
filename = args[1]

data <- read.table(filename)
a <- data$V2[1]
b <- data$V2[2]
c <- data$V2[3]
d <- data$V2[4]

results <- fisher.test(rbind(c(a,b),c(c,d)))

saveRDS(results, paste(filename, ".results.RDS", sep = ""))

sink(paste(filename, ".results", sep = ""))
print(results)
sink()
