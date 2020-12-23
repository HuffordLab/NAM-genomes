#!/usr/bin/env Rscript
library(DESeq2)
library(gplots)
library(calibrate)
library(RColorBrewer)
library(ape)

countdata <- read.table("NAMLINE_counts_cleaned-header.txt", sep="\t", header=TRUE, row.names=1)
coldata <- read.table("NAMLINE_metadata", sep="\t", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
dds = DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds = DESeq(dds)

png("NAMLINE_qc-dispersions.png", 1500, 1500, pointsize=20)
plotDispEsts(dds, main="NAMLINE Dispersion plot")
dev.off()

rld <- rlogTransformation(dds)

#hist(assay(rld))

condition <- coldata[,1]
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
sampleDists <- as.matrix(dist(t(assay(rld))))

png("NAMLINE_qc-heatmap-samples.png", w=1000, h=1000, pointsize=18)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"), dendrogram = c("row"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(14, 14), main="NAMLINE Sample Distance Matrix")
dev.off()

png("NAMLINE_qc-dendrogram-samples.png", w=1000, h=1000, pointsize=18)
sampledist <- dist(t(assay(rld)))
hc <- hclust(sampledist)
clus4 = cutree(hc, length(levels(condition))-2)
plot(as.phylo(hc), tip.color = mycols[clus4])
title("NAMLINE dendrogram")
dev.off()

rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="NAMLINE PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
#  rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
#  pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
#  terldt = list(levels(fac)), rep = FALSE)))
}

png("NAMLINE_qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition" )
dev.off()

res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"

#write.csv(resdata, file="NAMLINE-diffexpr-results.csv")
#attr(res, "filterThreshold")
#plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}

png("NAMLINE_diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="NAMLINE MA Plot")
dev.off()

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="NAMLINE Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

png("NAMLINE_diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()

