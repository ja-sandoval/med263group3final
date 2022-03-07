#BiocManager::install("tximportData")
#BiocManager::install("readr")
library(tximportData)
library(tximport)
library(readr)

### READ FILES
setwd('/Users/karlagodinez/Downloads/med263p')
# Experimental sample setup

samples <- read.table("samples_reduced.txt", header = TRUE)
samples

# ACCESS FILES
# files is a vector with the list of X RSEM output files
files <- file.path("rsem", paste0(samples$sample, ".genes.results"))
names(files) <- samples$sample


# Import files
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
names(txi)
head(txi$counts)


### DIFFERENTIALLY EXPRESSED GENES ANALYSIS
#https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#DESeq2
library(DESeq2)

# Define ctrl vs condition
sampleTable <- data.frame(condition = samples$pop,gender = samples$gender,replicate = samples$replicate)
sampleTable$condition <- relevel(factor(sampleTable$condition), ref = "healthy")
rownames(sampleTable) <- colnames(txi$counts)

# DEG
txi$length[txi$length <= 0] <- 1 # Non zero lengths
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

## Plots - Quality check
# Heatmap
library("pheatmap")
library("RColorBrewer")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd), vsd$condition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
         
# PCA
library(ggplot2)
plotPCA(vsd, intgroup = c('condition'))

pcaData <- plotPCA(vsd, intgroup='condition', returnData = T)
pcaData['gender'] <- vsd$gender
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + guides(color = guide_legend(order=2), shape = guide_legend(order=1)) + ggtitle("PCA NAFLD and NASH patients compared to healthy")
ggplot(pcaData, aes(x = PC1, y = PC2, shape = condition, color = gender)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + guides(color = guide_legend(order=2), shape = guide_legend(order=1)) + ggtitle("PCA NAFLD and NASH patients compared to healthy and gender")

# MA-plot
#BiocManager::install("apeglm")
library("apeglm")
resultsNames(dds)
res_maplot <- lfcShrink(dds, coef="condition_NAFLD_vs_healthy", type="apeglm")
par(mfrow=c(1,2))
plotMA(dds, main='Fold change vs mean expression\nNAFLD vs healthy')
plotMA(res_maplot, ylim = c(-5, 5),main='Normalized Fold change vs mean expression\nNAFLD vs healthy')
#idx <- identify(res_maplot$baseMean, res_maplot$log2FoldChange)
#rownames(res_maplot)[idx]
#dev.off()


# Conduct analysis

resultsNames(dds)
#res <- results(dds, contrast=c('condition','NAFLD','healthy'))
#resOrdered <- res[order(res$pvalue),]
#resOrdered <- res[order(res$log2FoldChange),]
resOrdered <- res_maplot[order(res_maplot$log2FoldChange),]


### GENE ID NAMES TO GENE SYMBOLS
#https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
library('biomaRt')

# Look for gene symbols #GTF 38 -v104
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(resOrdered)

# Keep genes with annotation
g_list <- getBM(filters= "ensembl_gene_id_version", attributes= c("ensembl_gene_id_version","hgnc_symbol",'gene_biotype'),values=genes,mart= mart)
g_list_na <- g_list[g_list$hgnc_symbol != '',]
g_list_na <- g_list_na[!duplicated(g_list_na[,2]),]

resOrdered_labeled <- subset(resOrdered, rownames(resOrdered) %in% g_list_na[,1])
rownames(resOrdered_labeled) <-  g_list_na$hgnc_symbol[ match(rownames(resOrdered_labeled), g_list_na$ensembl_gene_id_version) ]

# Add annotations -- gene symbol and type
resOrdered_labeled['gene_biotype'] <- g_list_na$gene_biotype[ match(rownames(resOrdered_labeled), g_list_na$ensembl_gene_id_version) ]
head(resOrdered_labeled)

### Significant Genes - Fold Change
sig_thresh = 2
resSig <-subset(resOrdered_labeled, padj < 0.05) # Significant
resUp <- subset(resOrdered_labeled, padj < 0.05 & log2FoldChange >= sig_thresh) # Significant & > threshold
resDown <- subset(resOrdered_labeled, padj < 0.05 & log2FoldChange <= -sig_thresh) # Significant & < -threshold

plot(x=resOrdered_labeled$log2FoldChange,y=resOrdered_labeled$padj)




