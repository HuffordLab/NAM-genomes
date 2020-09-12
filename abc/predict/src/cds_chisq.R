#NEW
#snp total 27256799
#number of cds snps 1322030
#sv total 98305
#sv total 37154

#snp total 27256799
#number of cds snps 1381301
#sv total 291539
#sv total 37154


snp_total <- 27256799
snp_cds <- 1322030
sv_total <- 98305
sv_cds <- 8686

count_mat <- matrix(
  c(snp_total - snp_cds, snp_cds,
    sv_total-sv_cds, sv_cds), ncol = 2,
  byrow = TRUE, dimnames = list(c("snp", "sv"), c("not cds", "cds"))
)

print(count_mat)

fisher.test(count_mat)

gs <- sum(count_mat)
rs <- rowSums(count_mat)
cs <- colSums(count_mat)
exp <- c(rs[1] * cs[1] / gs,  rs[1] * cs[2] / gs, rs[2] * cs[1] / gs, rs[2] * cs[2] / gs)
obs <- c(snp_total - snp_cds, snp_cds, sv_total-sv_cds, sv_cds )

print("((Obs - Exp)^2)/ Exp")
((obs - exp) ^ 2 ) / exp
obs
exp
pchisq(sum(((obs - exp) ^ 2 ) / exp), df = 1, lower.tail = F, log.p = T)
chisq.test(count_mat)
