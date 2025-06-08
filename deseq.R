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
Counts3g <- Counts %>% select(-c('ES_non_K1','ES_non_K2','ES_non_K3'))

# change column names
coldata <- data.frame(condition= factor(c("C","C","C","E","E","E","V","V","V","K","K","K")),row.names = c("C2","C3","C4","ES1","ES2","ES3","ES_non_K1","ES_non_K2","ES_non_K3","K1","K2","K3"))
     #check if the names and orders match
     all(colnames(Counts) %in% rownames(coldata))
     all(colnames(Counts) == rownames(coldata))
     
#For 3 group
coldata3g <- data.frame(condition= factor(c("C","C","C","E","E","E","K","K","K")),row.names = c("C2","C3","C4","ES1","ES2","ES3","K1","K2","K3"))
     #check if the names and orders match
     all(colnames(Counts3g) %in% rownames(coldata3g))
     all(colnames(Counts3g) == rownames(coldata3g))
     
#get DESeq matrix
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~ condition)

dds3g <- DESeqDataSetFromMatrix(countData = Counts3g, colData = coldata3g, design = ~ condition)

# filter out low expression gene
keep <- rowSums(counts(dds)) >= 12
dds <- dds[keep,]

keep <- rowSums(counts(dds3g)) >= 12
dds3g <- dds3g[keep,]

# set reference
dds$condition <- relevel(dds$condition, ref = "C")

dds3g$condition <- relevel(dds3g$condition, ref = "C")

#run DESeq
dds <- DESeq(dds, betaPrior = FALSE)

dds3g <- DESeq(dds3g, betaPrior = FALSE)

#quality control
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

#for 3 groups
vsdata3g <- vst(dds3g, blind = FALSE)
pcaData <- plotPCA(vsdata3g, intgroup=c("condition"), returnData=TRUE)
ggplot(pcaData, aes(PC1, PC2, color=condition), label = TRUE, label.size = 3) +
  geom_point(size=9) +
  coord_fixed()

# compare between samples, 2 at a time
res<- results(dds, alpha = 0.05, contrast = c("condition", "K", "C"))
summary(res)
plotMA(res)



res1<- results(dds, alpha = 0.05, contrast = c("condition", "E", "K"))
summary(res1)
plotMA(res1)

res2<- results(dds, alpha = 0.05, contrast = c("condition", "V", "c"))
summary(res2)
plotMA(res2)

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

sigs.df.raw <- as.data.frame(sigs)

sigsKvE.df <- as.data.frame(sigs1)

# filter the gene set to select the most differently expressed gene to the heatmap
sigs.df <- sigs.df.raw[(sigs.df.raw$baseMean > 100) & (abs(sigs.df.raw$log2FoldChange) > 1.6),]


sigs.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL", column = "SYMBOL")

# Filter out N/A symbol
sigs.df <- sigs.df[!is.na(sigs.df$symbol),]

dds$symbol <- mapIds(org.Mm.eg.db, keys = rownames(dds), keytype = "ENSEMBL", column = "SYMBOL")
#3groups heatmap
dds3g$symbol <- mapIds(org.Mm.eg.db, keys = rownames(dds3g), keytype = "ENSEMBL", column = "SYMBOL")

#heatmap, see arguments at https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")


library(ComplexHeatmap)

my_plot <- function(dds_in_function) {

  mat <- counts(dds_in_function, normalized = T)[rownames(sigs.df),]
  mat.z <- t(apply(mat, 1, scale))
  colnames(mat.z) <- rownames(coldata)
  sigs.df_filtered <- sigs.df[which(rownames(sigs.df) %in% matzgenes),]
  labels <- sigs.df$symbol
  
  hm.sigs <- Heatmap(mat.z, column_order = c('C2',"C3",'C4','K1','K2','K3','ES1','ES2','ES3','ES_non_K1','ES_non_K2','ES_non_K3'), 
                        #row_order = order(),
                        cluster_rows = T,name = "Z-score",
                        column_names_side = "top", column_labels = colnames(mat.z), column_names_rot = 30,
                        #column_split = factor(rep(c("Ctrl", "ES","Debr","ES_Debr"), 3)),
                        width = ncol(mat.z)*unit(10, "mm"),
                        height = nrow(mat.z)*unit(4, "mm")
  ) +
    rowAnnotation(labels = anno_text(labels, which = "row"),
                 width = max(grobWidth(textGrob(labels))
                  ))
  print('Modified hm.sigs')
  draw(hm.sigs, heatmap_legend_side = "left", gap = unit(0.1, "cm"))

}
my_plot(dds)

#to plot selective genes
my_plot_xjj <- function(dds_xjj) {
  
  matzgenes <- c('ENSMUSG00000038539','ENSMUSG00000064356','ENSMUSG00000069806','ENSMUSG00000004044',
                 'ENSMUSG00000000183','ENSMUSG00000025498',
                 'ENSMUSG00000036466','ENSMUSG00000050578','ENSMUSG00000005800','ENSMUSG00000032741',
                 'ENSMUSG00000032495','ENSMUSG00000039254','ENSMUSG00000044709','ENSMUSG00000036840',
                 'ENSMUSG00000029162','ENSMUSG00000075028','ENSMUSG00000039021','ENSMUSG00000020593','ENSMUSG00000033009',
                 'ENSMUSG00000020783','ENSMUSG00000041528','ENSMUSG00000028621','ENSMUSG00000088088','ENSMUSG00000040412',
                 'ENSMUSG00000065947','ENSMUSG00000056328','ENSMUSG00000044033',
                 'ENSMUSG00000041986','ENSMUSG00000029683')
  #labels <- c('Atf5','ATP8','Cacng7','Cavin1','Fgf6','Irf7','Megf11','Mmp13','Mmp8','Tpcn1','Lrrc2','Pomt1','Gemin7','Siah1a',
  #            'Khk','Prdm11','Ttc16','Lpin1','Ogfod1','Ncbp3','Rnf123','Cyb5rl','Rmrp','Elapor1','ND4L','Myh1','Ccdc141',
  #            'Elmod1','Lmod2') 
  
  sigsVvC.df1_filtered <- sigsVvC.df1[which(rownames(sigsVvC.df1) %in% matzgenes),]
  
  
  mat.xjj <- counts(dds_xjj, normalized = T)[rownames(sigsVvC.df1_filtered),]
  
  
  mat.z.xjj <- t(apply(mat.xjj, 1, scale))
  colnames(mat.z.xjj) <- rownames(coldata)
  
  
  
  
  mat.z.filtered <- mat.z.xjj
  

  
  # mat.z.filtered <- mat.z
  # labels <- sigsVvC.df1$symbol

  labels <- sigsVvC.df1_filtered$symbol
  hm.sigsVvC <- Heatmap(mat.z.filtered, column_order = c('C2',"C3",'C4','K1','K2','K3','ES1','ES2','ES3','ES_non_K1','ES_non_K2','ES_non_K3'), 
                        #row_order = order(),
                        cluster_rows = T, column_labels = colnames(mat.z.filtered),name = "Z-score",
                        column_names_side = "top",column_names_rot = 45,
                        #column_split = factor(rep(c("Ctrl", "ES","Debr","ES_Debr"), 3)),
                        width = ncol(mat.z.filtered)*unit(10, "mm"),
                        height = nrow(mat.z.filtered)*unit(3, "mm")
  ) +
    rowAnnotation(labels = anno_text(labels, which = "row"),
                  width = max(grobWidth(textGrob(labels))
                  ))
  
  draw(hm.sigsVvC, heatmap_legend_side = "left", gap = unit(0.1, "cm"))
  
}


my_plot_xjj(dds)

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

keyvals <- ifelse(
  res1$pvalue < 0.05 & res1$log2FoldChange > 1, 'red',
  ifelse(res1$pvalue < 0.05 & res1$log2FoldChange < -1, 'royalblue2',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'grey'] <- 'Non-sig'
names(keyvals)[keyvals == 'royalblue2'] <- 'Down'

EnhancedVolcano(res1,
                lab = res1$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('Map3k9',
                              'Akt1','Hspg2','Mbp','DUSP1','Slc6a12','Slc39a12'),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = "ES vs Sham",
                pCutoff = 0.05,
                pointSize = 2.0,
                labSize = 5.0,
                # col=c('grey86','deepskyblue2', 'chartreuse3', 'deeppink1'),
                colCustom = keyvals,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                legendLabels=c('Not sig.','Log (base 2) FC','adj-p',
                               'adj-p & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)

res$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res), keytype = "ENSEMBL", column = "SYMBOL")
res <- res[!is.na(res$symbol),]

keyvals1 <- ifelse(
  res$pvalue < 0.05 & res$log2FoldChange > 1, 'red',
  ifelse(res$pvalue < 0.05 & res$log2FoldChange < -1, 'royalblue2',
         'grey'))
keyvals1[is.na(keyvals1)] <- 'grey'
names(keyvals1)[keyvals1 == 'red'] <- 'Up'
names(keyvals1)[keyvals1 == 'grey'] <- 'Non-sig'
names(keyvals1)[keyvals1 == 'royalblue2'] <- 'Down'

EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('Atp8a2','Atf5',
                              'Akt1','Gpr37','Mbp','Hspg2'),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = "Control vs Keratectomy",
                pCutoff = 10e-2,
                pointSize = 2.0,
                labSize = 5.0,
                colCustom = keyvals1,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                legendLabels=c('Not sig.','Log (base 2) FC','adj-p',
                               'adj-p & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)

res2$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res2), keytype = "ENSEMBL", column = "SYMBOL")
res2 <- res2[!is.na(res2$symbol),]

keyvals2 <- ifelse(
  res2$pvalue < 0.05 & res2$log2FoldChange > 1, 'red',
  ifelse(res2$pvalue < 0.05 & res2$log2FoldChange < -1, 'royalblue2',
         'grey'))
keyvals2[is.na(keyvals2)] <- 'grey'
names(keyvals2)[keyvals2 == 'red'] <- 'Up'
names(keyvals2)[keyvals2 == 'grey'] <- 'Non-sig'
names(keyvals2)[keyvals2 == 'royalblue2'] <- 'Down'

EnhancedVolcano(res2,
                lab = res2$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('Mmp8',
                              'Mmp13','Megf11','Fgf6','Atf5','Cavin1','Tpcn1'),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = "ES only vs Control",
                pCutoff = 0.05,
                pointSize = 2.0,
                labSize = 5.0,
                colCustom = keyvals2,
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
david_KvE <- read.delim('C:/Users/yangm/R-Local/ES cornea/KvE_down.txt')

sig_filter <- function(df){
  # filter p.adj < 0.05
  sig_df <- df[df$Benjamini < 0.05,]
  # make a dot plot
  sig_df$Pathway1 <- sapply(strsplit(sig_df$Pathway,':'), tail, 1)
  
  return (sig_df)
}

david_KvE_sig <- sig_filter(david_KvE)


# c('Acute myeloid leukemia', 'Relaxin signaling pathway', 'Focal adhesion', 'Proteoglycans in cancer', 'Chemokine signaling pathway', 'Ras signaling pathway', 'PI3K-Akt signaling pathway', 'Disordered')
# c('Disordered','PI3K-Akt signaling pathway','Ras signaling pathway','Chemokine signaling pathway','Proteoglycans in cancer','Focal adhesion', 'Relaxin signaling pathway', 'Acute myeloid leukemia')
# david_KvE_sig$Pathway2 <- factor(david_KvE_sig$Pathway1, levels=c('PI3K-Akt signaling pathway','Ras signaling pathway','Chemokine signaling pathway','Proteoglycans in cancer','Focal adhesion', 'Relaxin signaling pathway', 'Acute myeloid leukemia'))

david_KvE_sig_KEGG <- david_KvE_sig %>% dplyr::filter(Category=='KEGG_PATHWAY') %>% arrange(Fold.Enrichment)

bubble_plot_pathway1 <- function(df) {
  df$Pathway1 <- factor(df$Pathway1, levels = df$Pathway1[order(df$Fold.Enrichment)])
  
  ggplot(df, aes(y = Pathway1, x = Fold.Enrichment, size=Count, color = Benjamini)) + 
    geom_point() +
    scale_size(range=c(10,20)) +
    scale_color_gradient(low = "blue1", high = "red1") +
    labs(color = "p-adj.") +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text = element_text(size = 25)) + 
    theme_classic(base_size = 20)
  
}

# bubble_plot(david_KvE_sig)
bubble_plot_pathway1(david_KvE_sig_KEGG)



#upregulated pathways
david_up <- read.delim('C:/Users/yangm/R-Local/ES cornea/KvE_up_modified.txt')
david_up_df <- david_up[david_up$PValue < 0.05,]
david_up_KEGG <- david_up_df %>% dplyr::filter(Category=='KEGG_PATHWAY') %>% arrange(Fold.Enrichment)
david_up_GOTERM <- david_up_df[startsWith(david_up_df$Category, "GOTERM"),] %>% arrange(Fold.Enrichment)

library(tools)

bubble_plot_term <- function(df) {
  df$Term <- toTitleCase(df$Term)
  df$Term <- factor(df$Term, levels = df$Term[order(df$Fold.Enrichment)])
  
  
  ggplot(df, aes(y = Term, x = Fold.Enrichment, size=Count, color = Benjamini)) + 
    geom_point() +
    scale_size(range=c(10,20)) +
    scale_color_gradient(low = "blue1", high = "red1") +
    labs(color = "p-adj.") +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text = element_text(size = 25)) + 
    theme_classic(base_size = 20)
  
}

bubble_plot_term(david_up_GOTERM)

#Cnetplot

library(clusterProfiler)
library(enrichplot)


