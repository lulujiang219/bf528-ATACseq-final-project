library("readr")
library("tidyr")
library(tidyverse)
library(DESeq2)

setwd("/projectnb/bf528/students/luluj/project_3_individual/rnaseq")


#DEseq
# Read the count matrix from output.txt
counts<- read.table("results/output.csv", sep = ",", header = TRUE, row.names = 1)

# Remove genes with zero counts across all samples
counts <- counts[rowSums(counts) > 0, ]

# Create the coldata data frame
coldata <- data.frame(samples = colnames(counts), 
                      condition = c(rep("CTL", 3), rep("KO", 3)),
                      row.names = "samples")
coldata$condition <- as.factor(coldata$condition)

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# Run the DESeq analysis
dds <- DESeq(dds)

# Extract the results comparing KO vs CTL
res <- results(dds, contrast = c("condition", "KO", "CTL"))

# Read the id2gene mapping file
id2gene <- read_delim("results/id2gene.txt", col_names = c("geneids", "genenames"))

# Combine the results with gene names, arrange by adjusted p-value, and select relevant columns
res_table <- res %>%
  as_tibble(rownames = "geneids") %>%
  left_join(id2gene, by = "geneids") %>%
  arrange(padj) %>%
  select(geneids, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, genenames)


# Count the number of differentially expressed genes at FDR threshold of 0.01
num_deg <- sum(res_table$padj < 0.01, na.rm = TRUE)

# Print the number of differentially expressed genes
cat("There are", num_deg, "differentially expressed genes at an FDR threshold of 0.01.")




#Histogram
library(ggplot2)

# Filter for significantly differentially expressed genes
sig_genes <- res_table %>% 
  filter(padj < 0.01)

# Histogram of log2 fold changes for significant genes
ggplot(sig_genes, aes(x = log2FoldChange)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  ggtitle("Histogram of log2 Fold Changes for Differentially Expressed Genes") +
  xlab("log2 Fold Change") +
  ylab("Frequency") +
  theme_minimal()




#volcano plot
library(ggplot2)
library(dplyr)

# Calculate significance: padj < 0.05
res_table <- res_table %>%
  mutate(Significant = padj < 0.05,
         Direction = if_else(log2FoldChange > 0, "Up", "Down"))

# Identify the top ten most significant genes
top_genes <- res_table %>%
  filter(Significant) %>%
  arrange(padj) %>%
  slice_head(n = 10)

# Create the volcano plot
volcano_plot <- ggplot(res_table, aes(x = log2FoldChange, y = -log10(pvalue), label = genenames)) +
  geom_point(aes(color = Significant), alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text(data = top_genes, aes(label = genenames), vjust = 2, hjust = 0.5, check_overlap = TRUE) +
  labs(title = "Volcano Plot of DE Genes",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal()

# Print the plot
print(volcano_plot)





#fgsea
library(fgsea)
library(readr)
library(dplyr)
library(qusage)
library(AnnotationDbi)
library(org.Hs.eg.db)



# loading gene sets
geneSets <- read.gmt("m2.all.v2023.2.Mm.symbols.gmt")
set.seed(42)  # for reproducibility

res_table <- res %>%
  as_tibble(rownames = "geneids") %>%
  left_join(id2gene, by = "geneids") %>%
  arrange(padj)

# Filtering DE genes based on a significance threshold (e.g., padj < 0.05)
significant_genes <- res_table %>% 
  filter(padj < 0.05) %>%
  pull(genenames)  # Extract gene names

# Check the top entries
head(significant_genes)


# Load enrichR if not already installed
if (!require("enrichR")) {
  install.packages("enrichR")
  library(enrichR)
}


databases <- enrichR::listEnrichrDbs()
mouse_databases <- databases[grep("Mouse", databases, ignore.case = TRUE)]
selected_databases <- mouse_databases$libraryName


enrich_res <- enrichR::enrichr(genes = significant_genes, databases = selected_databases)
enrichment_df <- enrich_res[[1]]
enrichment_df <- enrichment_df[order(enrichment_df$Adjusted.P.value), ]

# Visualize the top 10 results using ggplot2
library(ggplot2)
ggplot(head(enrichment_df, 10), aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_col(fill = 'steelblue') +
  labs(title = "Top Enriched Pathways", x = "Pathway", y = "-log10 Adjusted P-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_genes <- res_table %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(genenames, log2FoldChange) %>%
  deframe()  # This converts it to a named vector suitable for fgsea

set.seed(42)  # Ensuring reproducibility
all_genes_jittered <- all_genes + runif(length(all_genes), min = -1e-4, max = 1e-4)

# Run fgseaMultilevel with jittered stats
fgsea_results <- fgseaMultilevel(
  pathways = geneSets,
  stats = all_genes_jittered,
  minSize = 15,
  maxSize = 500
)
fgsea_results$abs_NES <- abs(fgsea_results$NES)
fgsea_results <- fgsea_results[order(fgsea_results$abs_NES, decreasing = TRUE),]
ggplot(head(fgsea_results, 20), aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_col() +
  coord_flip() +  # Horizontal bars
  scale_fill_manual(values = c("blue", "red"), 
                    labels = c("Downregulated", "Upregulated"), 
                    name = "Regulation") +
  labs(title = "Top FGSEA Results", x = "Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        legend.position = "top")