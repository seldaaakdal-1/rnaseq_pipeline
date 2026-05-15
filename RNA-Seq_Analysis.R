# ===============================
# RNA-seq Differential Expression Analysis Pipeline
# DESeq2 + Visualization + GO/KEGG Enrichment
# ===============================


# -------------------------------
# 1. Load required packages
# -------------------------------

# Load necessary Bioconductor and CRAN packages for RNA-seq analysis, differential expression, visualization, and functional enrichment.

# Bioconductor packages are used for RNA-seq and functional analysis
# Install if not already installed (commented for reproducibility)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install(c("DESeq2", "org.Sc.sgd.db", "clusterProfiler", "pathview", "EnhancedVolcano"))

# CRAN packages
# install.packages(c("tidyverse", "pheatmap"))


library(DESeq2)
library(org.Sc.sgd.db)
library(clusterProfiler)
library(pathview)
library(EnhancedVolcano)
library(pheatmap)
library(tidyverse)


# -------------------------------
# 2. Load count data
# -------------------------------

# Import raw gene count matrix generated from featureCounts.
count_data <- read.delim("gene_count_cutted.txt", sep="\t", row.names=1)

# -------------------------------
# 3. Define sample metadata
# -------------------------------

# Create experimental design table (condition: knockout vs wild type).
sample_info <- data.frame(
  condition = c("knockout", "knockout", "wild_type", "wild_type")
)

rownames(sample_info) <- colnames(count_data)

# Ensure sample names match count matrix columns.
all(rownames(sample_info) == colnames(count_data))


# -------------------------------
# 4. Create DESeq2 dataset
# -------------------------------

# Construct DESeq2 object using count data, metadata, and experimental design.

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ condition
)


# -------------------------------
# 5. Filtering low count genes
# -------------------------------

# Remove genes with very low read counts to improve statistical power.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# -------------------------------
# 6. Set reference level
# -------------------------------

# Define wild type as reference condition for differential expression analysis.
dds$condition <- relevel(dds$condition, ref = "wild_type")


# -------------------------------
# 7. Run DESeq2 analysis
# -------------------------------

# Perform normalization, dispersion estimation, and differential expression testing.
dds <- DESeq(dds)

# Extract results table (log2FC, p-values, adjusted p-values).
res <- results(dds)


# -------------------------------
# 8. PCA analysis
# -------------------------------

# Perform variance stabilizing transformation (VST).
vstdata <- vst(dds, blind = FALSE)

# Visualize sample clustering using PCA.
plotPCA(vstdata, intgroup = "condition")


# -------------------------------
# 9. MA plot
# -------------------------------

# Visualize log fold change versus mean expression for all genes.
plotMA(dds, alpha = 0.05)


# -------------------------------
# 10. Volcano plot preparation
# -------------------------------

res_df <- as.data.frame(res)

# Identify significantly differentially expressed genes.
res_df$significant <- ifelse(res_df$padj < 0.05, "yes", "no")

# Visualize significance (padj) vs fold change using volcano plot.
ggplot(res_df, aes(log(baseMean), log2FoldChange, color = significant)) +
  geom_point()


EnhancedVolcano(res_df,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(res_df))


# -------------------------------
# 11. Heatmap (Top genes)
# -------------------------------

vst_mat <- assay(vstdata)

# Select top significant genes based on p-value and fold change.
top_20_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(20) %>%
  rownames()

vst_top20 <- vst_mat[top_20_genes, ]
colnames(vst_top20) <- sample_info$condition

# Visualize expression patterns using scaled heatmap.
pheatmap(vst_top20, 
         scale = "row",                     
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Heatmap Plot")


# -------------------------------
# 12. GO enrichment analysis
# -------------------------------

# Identify enriched Gene Ontology terms (BP, MF, CC)
# Using significantly expressed genes.

res_df$significant <- ifelse(res_df$padj < 0.05, "yes", "no")
res_go <- res_df %>%
  filter(significant == "yes")

gene_list <- rownames(res_go) 
head(gene_list)

ego1 <- enrichGO(gene = gene_list,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF", 
                 ont = "BP")

barplot(ego1, title = "Biological Process")

ego2 <- enrichGO(gene = gene_list,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF",
                 ont = "MF")

barplot(ego2, title = "Molecular Function")

ego3 <- enrichGO(gene = gene_list,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF",
                 ont = "CC")

ego3_plot <- barplot(ego3, title = "Cellular Component") 
ego3_plot


# -------------------------------
# 13. KEGG pathway analysis
# -------------------------------

search <- search_kegg_organism("Saccharomyces cerevisiae") 

search 

# Perform pathway enrichment analysis using KEGG database.
kegg_enrich <- enrichKEGG(gene = gene_list,
                          organism = "sce",
                          keyType = "kegg")

head(kegg_enrich)

# Visualize significant pathways and map gene expression changes.
browseKEGG(kegg_enrich, "sce00190") 

gen_verisi <- res_go$log2FoldChange

names(gen_verisi) <- rownames(res_go)

sce00190<-pathview(gene.data = gen_verisi, 
                   species = "sce", 
                   pathway.id = "sce00190", 
                   gene.idtype = "orf")