library(tximportData)
# dir <- system.file("extdata", package = "tximportData")
# list.files(dir)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Step 1: Prepare Files

library(GenomicFeatures)

# Replace with your local GTF file path (mouse genome, e.g., GRCm39)
txdb <- makeTxDbFromGFF("Mus_musculus.GRCm39.109.gtf.gz", format = "gtf")

# Extract transcript-to-gene mapping
tx2gene <- select(txdb, keys = keys(txdb, keytype = "TXNAME"),
                  columns = c("GENEID"), keytype = "TXNAME")



library(tximport)
library(readr)
library(DESeq2)

# Step 2: Import Quantifications
samples <- read.csv("metadata.csv")
files <- setNames(samples$filepath, samples$sample)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Step 3: Normalize Using DESeq2
# Create a simple sample metadata table for DESeq2
sampleTable <- data.frame(condition = factor(c("C1","C2","C3","C4","K1","K2","K3","K4", "ES1","ES2","ES3","ES4","ES_non_K1", "ES_non_K2", "ES_non_K3")))
rownames(sampleTable) <- names(files)

# Build DESeq2 object
# dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~1)

# Perform normalization
# dds <- estimateSizeFactors(dds)
# normalized_counts <- counts(dds, normalized = TRUE)

# Keep rowname

#write.csv(normalized_counts, file = "normalized_counts.csv")

# Output to CSV, convert rowname to column "gene_id"

# # Convert matrix to data frame and add gene IDs as a column
# normalized_df <- as.data.frame(normalized_counts)
# normalized_df$gene_id <- rownames(normalized_counts)
# 
# # Move gene_id to the first column
# normalized_df <- normalized_df[, c("gene_id", setdiff(names(normalized_df), "gene_id"))]
# 
# # Write to CSV
# write.csv(normalized_df, file = "normalized_counts.csv", row.names = FALSE)

# Get DEG
# Define another sample_list to replace sampleTable below
sample_list <- c("C1","C2","C3","C4", "ES1","ES2","ES3","ES4","K1","K2","K3","K4","ES_non_K1", "ES_non_K2", "ES_non_K3")

# Create condition vector
condition <- gsub("[0-9]", "", sample_list)

# Build colData
colData <- data.frame(
  row.names = sample_list,
  condition = factor(condition, levels = c("C", "ES","K","ES_non_K" ))  # C is reference
)
dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ condition)
dds <- DESeq(dds)

#quality control
library(ggplot2)

vsdata <- vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "condition", returnData = TRUE)

plotDispEsts(dds)

#plot PCA by hand to make it prettier
pcaData <- plotPCA(vsdata, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=9) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# Remove C1, ES4 & K4 from counts since they seem to be outlier from plotPCA
library('dplyr')
samples_to_remove <- c("C1", "ES4", "K4")

dds_filtered <- dds[, !(colnames(dds) %in% samples_to_remove) ]
# Suppose your design is ~ group
dds_filtered$condition <- droplevels(dds_filtered$condition)

# Re-run DESeq
dds_filtered <- DESeq(dds_filtered)
#Re-run PCA
vsdata <- vst(dds_filtered, blind = FALSE)
plotPCA(vsdata, intgroup = "condition", returnData = TRUE)

plotDispEsts(dds_filtered)

#plot PCA by hand to make it prettier
pcaData <- plotPCA(vsdata, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=9) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

res <- results(dds_filtered, contrast = c("condition", "C", "K"))
write.csv(as.data.frame(res), file = "DEG.csv")

# change Ensembl id to gene symbol
# # Install if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(biomaRt)

# Connect to Ensembl BioMart for mouse
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ensembl_ids <- rownames(res)

mapping <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
# Prepare DESeq2 results for merging
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- sub("\\..*", "", rownames(res))  # cleaned Ensembl IDs

# Merge with gene symbol annotations
res_annotated <- merge(res_df, mapping, by = "ensembl_gene_id", all.x = TRUE)

# Optional: reorder columns
res_annotated <- res_annotated[, c("ensembl_gene_id", "mgi_symbol", setdiff(colnames(res_annotated), c("ensembl_gene_id", "mgi_symbol")))]

# Save to file
write.csv(res_annotated, file = "DEG_with_symbols_mouse.csv", row.names = FALSE)

# take out the significant gene
sigs <- na.omit(res_annotated)
sigs <- sigs[sigs$pvalue <0.05,]
sigs.df.raw <- as.data.frame(sigs)

sigs_adj <- sigs[sigs$padj <0.05,]
sigs_adj.df <- as.data.frame(sigs_adj)

# Volcano plot
library("ggrepel")
library("EnhancedVolcano")

keyvals <- ifelse(
  res_annotated$pvalue < 0.05 & res_annotated$log2FoldChange > 1, 'red',
  ifelse(res_annotated$pvalue < 0.05 & res_annotated$log2FoldChange < -1, 'royalblue2',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'grey'] <- 'Non-sig'
names(keyvals)[keyvals == 'royalblue2'] <- 'Down'

EnhancedVolcano(res_annotated,
                lab = res_annotated$mgi_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-15, 15), 
                selectLab = c( "Ninj2","Pomc","Otx2","Neu4", "Npb", "Ntsr2", "Akt1","Nnat","Il22ra2", "Nlrp4f", "Tlr8", "Il31ra"),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = "Control vs Injury",
                pCutoff = 0.05,
                pointSize = 2.0,
                labSize = 5.0,
                colCustom = keyvals,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                #legendLabels=c('Not sig.','Log (base 2) FC','p',
                #               'p & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)

#make heatmap
library(ComplexHeatmap)
# install.packages("circlize")
# install.packages("RColorBrewer")

library(circlize)
library(RColorBrewer)

final_sample_list <- c("C2","C3","C4", "ES1","ES2","ES3","K1","K2","K3","ES_non_K1", "ES_non_K2", "ES_non_K3")
condition <- gsub("[0-9]", "", final_sample_list)

coldata <- data.frame(
  row.names = final_sample_list,
  condition = factor(condition, levels = c("C", "ES","K","ES_non_K"))  # C is reference
)
sigs_adj_filtered <- sigs_adj.df[sigs_adj.df$log2FoldChange < -2.5 | sigs_adj.df$log2FoldChange > 2, ]
# Define genes to exclude explicitly
excluded_genes <- c("Dlk1", "Pitx1", "Cga", "Gh", "Prl", "Lhb", "Pomc", "Ltf","Mpo","Sh2d1a","Pramel56")

# Remove genes starting with "Gm" or ending with "Rik" in the mgi_symbol column
sigs_adj_filtered <- sigs_adj_filtered[ 
  !grepl("^Gm", sigs_adj_filtered$mgi_symbol) & 
    !grepl("Rik$", sigs_adj_filtered$mgi_symbol) &
    !(sigs_adj_filtered$mgi_symbol %in% excluded_genes),
]
rownames(sigs_adj_filtered) <- sigs_adj_filtered$ensembl_gene_id

my_plot <- function(dds_in_function) {
  
  mat <- counts(dds_in_function, normalized = T)[rownames(sigs_adj_filtered),]
  mat.z <- t(apply(mat, 1, scale))
  colnames(mat.z) <- rownames(coldata)
  labels <- sigs_adj_filtered$mgi_symbol
  
  hm.sigs <- Heatmap(mat.z, column_order = c("C2","C3","C4","K1","K2","K3", "ES1","ES2","ES3","ES_non_K1", "ES_non_K2", "ES_non_K3"), 
                     #row_order = order(),
                     cluster_rows = T,name = "Z-score",
                     column_names_side = "top", column_labels = colnames(mat.z), column_names_rot = 30,
                     width = ncol(mat.z)*unit(10, "mm"),
                     height = nrow(mat.z)*unit(4, "mm")
  ) +
    rowAnnotation(labels = anno_text(labels, which = "row"),
                  width = max(grobWidth(textGrob(labels))
                  ))
  print('Modified hm.sigs')
  draw(hm.sigs, heatmap_legend_side = "left", gap = unit(0.1, "cm"))
  
}
my_plot(dds_filtered)

#Now compare ES to K

res_ES <- results(dds_filtered, contrast = c("condition", "K", "ES"))
write.csv(as.data.frame(res_ES), file = "DEG_ES.csv")

# change Ensembl id to gene symbol
# # Install if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")

#library(biomaRt)

# Connect to Ensembl BioMart for mouse
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ensembl_ids <- rownames(res_ES)

mapping <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
# Prepare DESeq2 results for merging
res_ES_df <- as.data.frame(res_ES)
res_ES_df$ensembl_gene_id <- sub("\\..*", "", rownames(res))  # cleaned Ensembl IDs

# Merge with gene symbol annotations
res_ES_annotated <- merge(res_ES_df, mapping, by = "ensembl_gene_id", all.x = TRUE)

# Optional: reorder columns
res_ES_annotated <- res_ES_annotated[, c("ensembl_gene_id", "mgi_symbol", setdiff(colnames(res_annotated), c("ensembl_gene_id", "mgi_symbol")))]

# Save to file
#write.csv(res_ES_annotated, file = "DEG_ES_symbol.csv", row.names = FALSE)

# take out the significant gene
sigs_ES <- na.omit(res_ES_annotated)
sigs_ES <- sigs_ES[sigs_ES$pvalue <0.05,]
sigs_ES.df <- as.data.frame(sigs_ES)

sigs_ES_adj <- sigs_ES[sigs_ES$padj <0.05,]
sigs_ES_adj.df <- as.data.frame(sigs_ES_adj)

# Volcano plot
library("ggrepel")
library("EnhancedVolcano")

keyvals <- ifelse(
  res_ES_annotated$pvalue < 0.05 & res_ES_annotated$log2FoldChange > 1, 'red',
  ifelse(res_ES_annotated$pvalue < 0.05 & res_ES_annotated$log2FoldChange < -1, 'royalblue2',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'grey'] <- 'Non-sig'
names(keyvals)[keyvals == 'royalblue2'] <- 'Down'

EnhancedVolcano(res_ES_annotated,
                lab = res_ES_annotated$mgi_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-15, 15), 
                selectLab = c( "Cacna1s","Pomc", "Wnt3", "Wnt4", "Hif3a", "Muc6", "Bmp7", "Wnt11", "Sox9", "Bmp4"),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = "ES vs Injury",
                pCutoff = 0.05,
                pointSize = 2.0,
                labSize = 5.0,
                colCustom = keyvals,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                #legendLabels=c('Not sig.','Log (base 2) FC','p',
                #               'p & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)

#Now plot ES_non_K vs Control
res_NK <- results(dds_filtered, contrast = c("condition", "C", "ES_non_K"))

# change Ensembl id to gene symbol
# # Install if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")

#library(biomaRt)

# Connect to Ensembl BioMart for mouse
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ensembl_ids <- rownames(res_NK)

mapping <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
# Prepare DESeq2 results for merging
res_NK_df <- as.data.frame(res_NK)
res_NK_df$ensembl_gene_id <- sub("\\..*", "", rownames(res))  # cleaned Ensembl IDs

# Merge with gene symbol annotations
res_NK_annotated <- merge(res_NK_df, mapping, by = "ensembl_gene_id", all.x = TRUE)

# Optional: reorder columns
res_NK_annotated <- res_NK_annotated[, c("ensembl_gene_id", "mgi_symbol", setdiff(colnames(res_NK_annotated), c("ensembl_gene_id", "mgi_symbol")))]

# Save to file
#write.csv(res_ES_annotated, file = "DEG_ES_symbol.csv", row.names = FALSE)

# take out the significant gene
sigs_NK <- na.omit(res_NK_annotated)
sigs_NK <- sigs_NK[sigs_NK$pvalue <0.05,]

sigs_NK_adj <- sigs_NK[sigs_NK$padj <0.05,]
sigs_NK_adj.df <- as.data.frame(sigs_NK_adj)

# Volcano plot
# library("ggrepel")
# library("EnhancedVolcano")

keyvals <- ifelse(
  res_NK_annotated$pvalue < 0.05 & res_NK_annotated$log2FoldChange > 1, 'red',
  ifelse(res_NK_annotated$pvalue < 0.05 & res_NK_annotated$log2FoldChange < -1, 'royalblue2',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'grey'] <- 'Non-sig'
names(keyvals)[keyvals == 'royalblue2'] <- 'Down'

EnhancedVolcano(res_NK_annotated,
                lab = res_NK_annotated$mgi_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-15, 15), 
                selectLab = c( "Cacna1s","Pomc", "Kcna7", "Tlr8", "Kcnj15", "Calcrl", "Calcoco2", "Ltb4r1", "Cx3cl1", "Kcng1", "Cacng3", "Kcnn4"),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = "Control vs ES without Injury",
                pCutoff = 0.05,
                pointSize = 2.0,
                labSize = 5.0,
                colCustom = keyvals,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                #legendLabels=c('Not sig.','Log (base 2) FC','p',
                #               'p & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)



# import From David
GO_BP <- read.delim('GO_BP.txt')
GO_BP$Term <- sub("GO:\\d+~", "", GO_BP$Term)

terms_to_include <- c(
  "calcium ion transport", "chemical synaptic transmission", 
  "transforming growth factor beta receptor signaling pathway",
  "cellular response to mechanical stimulus", "positive regulation of interleukin-6 production",
  "neuron projection development", "axon guidance",
  "potassium ion transmembrane transport", "insulin-like growth factor receptor signaling pathway",
  "positive regulation of neuron differentiation", "cell fate commitment", "synapse organization",
  "sodium ion transport", "canonical Wnt signaling pathway",
  "ERK1 and ERK2 cascade", "brain-derived neurotrophic factor receptor signaling pathway",
  "cellular oxidant detoxification", "response to retinoic acid", "axonogenesis",
  "macrophage colony-stimulating factor signaling pathway", "oligodendrocyte differentiation",
  "regulation of neuron projection development", "neuromuscular junction development", 
  "sodium ion transmembrane transport","synaptic membrane adhesion", "regulation of presynaptic cytosolic calcium ion concentration",
  "myelination in peripheral nervous system", "neurotransmitter transport",
  "neurofilament cytoskeleton organization"
  )

# GO_BP_filtered <- GO_BP[GO_BP$Benjamini < 0.01,]
# GO_BP_filtered <- GO_BP_filtered[GO_BP_filtered$Fold.Enrichment > 3, ]

GO_BP_filtered <- GO_BP_filtered[GO_BP_filtered$Term %in% terms_to_include, ]
GO_BP_filtered<- GO_BP_filtered[startsWith(GO_BP_filtered$Category, "GOTERM"),] %>% arrange(Fold.Enrichment)

library(tools)

bubble_plot_term <- function(df) {
  df$Term <- toTitleCase(df$Term)
  df$Term <- factor(df$Term, levels = df$Term[order(df$Fold.Enrichment)])
  
  
  ggplot(df, aes(y = Term, x = Fold.Enrichment, size=Count, color = Benjamini)) + 
    geom_point() +
    scale_size(range=c(3,8)) +
    scale_color_gradient(low = "blue1", high = "red1") +
    labs(color = "p-adj.") +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text = element_text(size = 25)) + 
    theme_classic(base_size = 20)
  
}

bubble_plot_term(GO_BP_filtered)

KEGG <- read.delim('KEGG.txt')

KEGG$Term <- sub("^[^:]+:", "", KEGG$Term)
KEGG_filtered <- KEGG[KEGG$Benjamini < 0.05,]
#GO_BP_filtered <- GO_BP_filtered %>% dplyr::filter(Category=='GOTERM') %>% arrange(Fold.Enrichment)
KEGG_filtered<- KEGG_filtered[startsWith(KEGG_filtered$Category, "KEGG_PATHWAY"),] %>% arrange(Fold.Enrichment)

library(tools)

bubble_plot_term <- function(df) {
  df$Term <- toTitleCase(df$Term)
  df$Term <- factor(df$Term, levels = df$Term[order(df$Fold.Enrichment)])
  
  
  ggplot(df, aes(y = Term, x = Fold.Enrichment, size=Count, color = Benjamini)) + 
    geom_point() +
    scale_size(range=c(4,12)) +
    scale_color_gradient(low = "blue1", high = "red1") +
    labs(color = "p-adj.") +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text = element_text(size = 25)) + 
    theme_classic(base_size = 20)
  
}

bubble_plot_term(KEGG_filtered)
