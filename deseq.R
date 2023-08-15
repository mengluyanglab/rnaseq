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
Counts <- read.csv("//nas.ltc.upinthecloud.info/homes/RNASeq/matrix/matrix.csv", header = TRUE, row.names = 1, sep = "\t")

# Remove C1 & K4 from counts since they seem to be outlier from plotPCA
library('dplyr')
#Counts <- Counts %>% select(-c('C1','K4'))

# change column names
coldata <- data.frame(condition= factor(c("C","C","C","C","E","E","E","E","V","V","V","K","K","K","K")),row.names = c("C1","C2","C3","C4","ES1","ES2","ES3","ES4","ES_non_K1","ES_non_K2","ES_non_K3","K1","K2","K3","K4"))
     #check if the names and orders match
     all(colnames(Counts) %in% rownames(coldata))
     all(colnames(Counts) == rownames(coldata))

#get DESeq matrix
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~ condition)

# filter out low expression gene
keep <- rowSums(counts(dds)) >= 15
dds <- dds[keep,]


# set reference
dds$condition <- relevel(dds$condition, ref = "C")

#run DESeq
dds <- DESeq(dds, betaPrior = FALSE)

#quality control
vsdata <- vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "condition", returnData = TRUE)

plotDispEsts(dds)

# compare between samples, 2 at a time
res<- results(dds, alpha = 0.05, contrast = c("condition", "K", "C"))
summary(res)
plotMA(res)

# draw a heatmap
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

dds$symbol <- mapIds(org.Mm.eg.db, keys = rownames(dds), keytype = "ENSEMBL", column = "SYMBOL")

#heatmap
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")


library(ComplexHeatmap)

mat <- counts(dds, normalized = T)[rownames(sigs.df),]


mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(coldata)


labels <- sigs.df$symbol


hm.sigs <- Heatmap(mat.z, column_order = c('C2',"C3",'C4','K1','K2','K3','ES1','ES2','ES3','ES4','ES_non_K1','ES_non_K2','ES_non_K3'), cluster_rows = T, column_labels = colnames(mat.z),name = "Z-score") +
  rowAnnotation(labels = anno_text(labels, which = "row"), 
                width = max(grobWidth(textGrob(labels))))
draw(hm.sigs, gap = unit(1, "cm"))

# add annotation "https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html"

sigsKvsE.df <- as.data.frame(sigsKvsE)
sigsKvsE.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sigsKvsE.df), keytype = "ENSEMBL", column = "SYMBOL")

matKvsE <- counts(dds, normalized = T)[rownames(sigsKvsE.df),]


matKvsE.z <- t(apply(matKvsE, 1, scale))
colnames(matKvsE.z) <- rownames(coldata)


labels <- sigsKvsE.df$symbol


hm.sigs <- Heatmap(matKvsE.z, column_order = c('C2',"C3",'C4','K1','K2','K3','ES1','ES2','ES3','ES4','ES_non_K1','ES_non_K2','ES_non_K3'), cluster_rows = T, column_labels = colnames(matKvsE.z),name = "Z-score") +
  rowAnnotation(labels = anno_text(labels, which = "row"), 
                width = max(grobWidth(textGrob(labels))))
draw(hm.sigs, gap = unit(1, "cm"))



