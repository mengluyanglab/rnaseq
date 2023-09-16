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
install.packages("dplyr")
library('dplyr')
Counts <- Counts %>% select(-c('C1','K4','ES4'))

# change column names
coldata <- data.frame(condition= factor(c("C","C","C","E","E","E","V","V","V","K","K","K")),row.names = c("C2","C3","C4","ES1","ES2","ES3","ES_non_K1","ES_non_K2","ES_non_K3","K1","K2","K3"))
     #check if the names and orders match
     all(colnames(Counts) %in% rownames(coldata))
     all(colnames(Counts) == rownames(coldata))

#get DESeq matrix
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~ condition)

# filter out low expression gene
keep <- rowSums(counts(dds)) >= 12
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

res1<- results(dds, alpha = 0.05, contrast = c("condition", "E", "K"))
summary(res1)
plotMA(res1)

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

sigsKvE.df <- as.data.frame(sigs1)

# filter the gene set to select the most differently expressed gene to the heatmap
sigs.df1 <- sigs.df[(sigs.df$baseMean > 100) & (abs(sigs.df$log2FoldChange) > 1.75),]


sigs.df1$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sigs.df1), keytype = "ENSEMBL", column = "SYMBOL")

# Filter out N/A symbol
sigs.df1 <- sigs.df1[!is.na(sigs.df1$symbol),]

dds$symbol <- mapIds(org.Mm.eg.db, keys = rownames(dds), keytype = "ENSEMBL", column = "SYMBOL")

#heatmap, see arguments at https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")


library(ComplexHeatmap)

mat <- counts(dds, normalized = T)[rownames(sigsKvE.df1),]


mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(coldata)


labels <- sigsKvE.df1$symbol


hm.sigsVvC <- Heatmap(matVvC.z, column_order = c('C2',"C3",'C4','ES_non_K1','ES_non_K2','ES_non_K3','K1','K2','K3','ES1','ES2','ES3'), 
                      cluster_rows = T, column_labels = colnames(matVvC.z),name = "Z-score",
                      column_names_side = "top",column_names_rot = 45,
                      #column_split = factor(rep(c("Ctrl", "ES","Debr","ES_Debr"), 3)),
                      width = ncol(matVvC.z)*unit(10, "mm"),
                      height = nrow(matVvC.z)*unit(4, "mm")
) +
  rowAnnotation(labels = anno_text(labels, which = "row"),
               width = max(grobWidth(textGrob(labels))
                ))

draw(hm.sigsVvC, gap = unit(0.1, "cm"))

#Making volcano plot using (https://github.com/kevinblighe/EnhancedVolcano)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
# this depend on ggplot2 and ggrepel, so install ggrepel, and read ggplot2
install.packages("ggrepel")

library("ggrepel")
library("EnhancedVolcano")

#Annotate gene symbols to res
res1$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res1), keytype = "ENSEMBL", column = "SYMBOL")
res1 <- res1[!is.na(res1$symbol),]


EnhancedVolcano(res1,
                lab = res1$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('Camk1','Map3k9','Il17rb',
                              'Akt1','Hspg2','Mbp','DUSP1','Slc6a12','Slc39a12'),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = "ES vs Sham",
                pCutoff = 10e-2,
                pointSize = 2.0,
                labSize = 5.0,
                col=c('black', 'black', 'black', 'red2'),
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                legendLabels=c('Not sig.','Log (base 2) FC','adj-p',
                               'adj-p & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)

# add annotation "https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html"

# get row na   .names(sigsKvE.df),file='sigsKvE-rownames.csv',quote=FALSE,row.names=FALSE, col.names=FALSE)

# GO 
BiocManager::install("clusterProfiler")

library(clusterProfiler)
library("AnnotationDbi")
library("org.Mm.eg.db")

if (!("msigdbr" %in% installed.packages())) {
  BiocManager::install("msigdbr", update = FALSE)
}

library(msigdbr)

# To use %>%
install.packages("magrittr")
library(magrittr)

msigdbr_species()
mm_hallmark_sets <- msigdbr(species = "Mus musculus", category = "H")
mm_C5_sets <- msigdbr(species = "Mus musculus", category = "C5")

head(mm_hallmark_sets)
keytypes(org.Mm.eg.db)

# Create a mapped data frame
dge_mapped_df <- data.frame(
    gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Mm.eg.db,
    keys = rownames(sigs.df),
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
dplyr::filter(!is.na(gene_symbol)) %>%
# Make an `Ensembl` column to store the rownames
tibble::rownames_to_column("Ensembl") %>% 
## Now let's join the rest of the expression data
dplyr::inner_join(tibble::rownames_to_column(sigs.df, var="Ensembl"), by = "Ensembl")

#named vector ranked based on the log2 fold change values
lfc_vector <- dge_mapped_df$log2FoldChange
names(lfc_vector) <- dge_mapped_df$gene_symbol

# sort log2foldchange in decending order
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# Set the seed so our results are reproducible:
set.seed(2023)

gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 20000, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_C5_sets,
    gs_name,
    gene_symbol
  )
)

head(gsea_results@result)
gsea_result_df <- gsea_results@result
gsea_result_sig <- gsea_result_df[gsea_result_df$p.adjust < 0.2,]

#Find upregulation gene list
genes_to_test <- rownames(sigs[sigs$log2FoldChange > 0,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)

#fit <- plot(barplot(GO_results, showCategory = 20))

# import From David
david_KvE <- read.delim('C:/Users/yangm/R-Local/rnaseq/KvE_down.txt')
# filter p.adj < 0.05
david_KvE_sig <- david_KvE[david_KvE$Benjamini < 0.05,]
# make a dot plot
david_KvE_sig$Pathway1 <- sapply(strsplit(david_KvE_sig$Pathway,':'), tail, 1)
# c('Acute myeloid leukemia', 'Relaxin signaling pathway', 'Focal adhesion', 'Proteoglycans in cancer', 'Chemokine signaling pathway', 'Ras signaling pathway', 'PI3K-Akt signaling pathway', 'Disordered')
# c('Disordered','PI3K-Akt signaling pathway','Ras signaling pathway','Chemokine signaling pathway','Proteoglycans in cancer','Focal adhesion', 'Relaxin signaling pathway', 'Acute myeloid leukemia')
david_KvE_sig$Pathway2 <- factor(david_KvE_sig$Pathway1, levels=c('Disordered','PI3K-Akt signaling pathway','Ras signaling pathway','Chemokine signaling pathway','Proteoglycans in cancer','Focal adhesion', 'Relaxin signaling pathway', 'Acute myeloid leukemia'))

ggplot(david_KvE_sig, aes(y = Pathway2, x = Fold.Enrichment, size=Count, color = Benjamini)) + 
  geom_point() +
  scale_size(range=c(6,10)) +
  labs(color = "p-adj.") +
  theme(axis.title.y=element_blank())


 
