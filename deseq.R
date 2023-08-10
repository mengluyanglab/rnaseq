# Install bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

# install DEseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

#dilate the first row if it is comments
Counts <- read.delim("//nas/homes/RNASeq/matrix/matrix.csv", header = TRUE, row.names = 1, sep = "\t")

# filter out low expression gene

filtered <- Counts[which(rowSums(Counts) >50),]

# change column name
filtered <- setNames(filtered, c("C1","C2","C3","C4","ES1","ES2","ES3","ES4","ESnonK1","ESnonK2","ESnonK3","K1","K2","K3","K4"))
     
# define which is control/treatment
condition <- factor(c("C","C","C","C","E","E","E","E","V","V","V","K","K","K","K"))

coldata <- data.frame(row.names = colnames(filtered),condition)

dds <- DESeqDataSetFromMatrix(countData = filtered, colData = coldata, design = ~condition)

dds <- DESeq(dds)

vsdata <- vst(dds, blind = FALSE)

#quality control
plotPCA(vsdata, intgroup = "condition")

plotDispEsts(dds)

# compare between samples, 2 at a time
res <- results(dds, contrast = c("condition","E","K"))

# take out the significant gene
sigs <- na.omit(res)
sigs <- sigs[sigs$padj <0.05,]

# convert Ensemble id to gene name
if (!require("BiocManager", quietly = TRUE))
  +   install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")

if (!require("BiocManager", quietly = TRUE))
  +   install.packages("BiocManager")
BiocManager::install("AnnotationDbi")

library("AnnotationDbi")
library("org.Mm.eg.db")

sigs.df <- as.data.frame(sigs)
sigs.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL", column = "SYMBOL")

filtered$symbol <- mapIds(org.Mm.eg.db, keys = rownames(filtered), keytype = "ENSEMBL", column = "SYMBOL")

#heatmap
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")


library(ComplexHeatmap)

mat <- counts(dds, normalized = T)[rownames(sigs.df),]


mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(coldata)


labels <- sigs.df$symbol


hm.sigs <- Heatmap(mat.z, column_order = rownames(coldata), cluster_rows = T, column_labels = colnames(mat.z),name = "Z-score") +
  rowAnnotation(labels = anno_text(labels, which = "row"), 
                width = max(grobWidth(textGrob(labels))))
draw(hm.sigs, gap = unit(1, "cm"))

# add annotation "https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html"





