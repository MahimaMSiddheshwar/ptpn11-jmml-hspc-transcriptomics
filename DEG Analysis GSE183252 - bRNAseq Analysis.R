# ==============================================================================
# JMML vs Normal HSPCs - GSE183252
# Bulk RNA-seq End-to-End Pipeline using DESeq2
# ==============================================================================

# ---- Setup & Libraries ----
library(GEOquery)
library(DESeq2)
library(apeglm)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(EnsDb.Hsapiens.v86)

# Avoid namespace clashes with ensembldb
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename

# ---- Settings ----
GEO_ID      <- "GSE183252"
OUTPUT_DIR  <- file.path("results", paste0(GEO_ID, "_DESeq2"))
QVAL_CUTOFF <- 0.05
FC_CUTOFF   <- 1.2
LFC_CUTOFF  <- log2(FC_CUTOFF)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- 1. Download Data ----
message("Downloading data from GEO...")
gse <- getGEO(GEO_ID, GSEMatrix = TRUE)[[1]]
supp <- getGEOSuppFiles(GEO_ID)

counts_file <- rownames(supp)[grepl("count", rownames(supp), ignore.case = TRUE)][1]
raw_counts <- read.table(counts_file, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)
raw_counts <- as.matrix(raw_counts)


# ---- 2. Sample Metadata Preparation ----
sample_info <- pData(gse)
sample_info$disease <- ifelse(grepl("JMML", sample_info$`disease state:ch1`),
                              "JMML", "Normal")

rownames(sample_info) <- sample_info$title
common <- intersect(colnames(raw_counts), rownames(sample_info))

# DESeq2 REQUIRES raw, un-normalized integer counts
counts <- round(raw_counts[, common])

coldata <- data.frame(
  disease   = factor(sample_info[common, "disease"], levels = c("Normal", "JMML")),
  row.names = common
)

print(table(coldata$disease))


# ---- 3. Pre-filtering Low Counts ----.
# Keep genes that have at least 10 reads in at least 15 samples.
keep_expr <- rowSums(counts >= 10) >= 15
counts_filt <- counts[keep_expr, ]


# ---- 4. DESeq2 Analysis Loop ----
message("Running DESeq2 differential expression pipeline...")

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                              colData = coldata,
                              design = ~ disease)

# Ensure "Normal" is treated as the reference control baseline
dds$disease <- relevel(dds$disease, ref = "Normal")

# Run the complete DESeq2 workflow (Normalization, Dispersions, Wald Test GLM)
dds <- DESeq(dds)

# Apply Modern Log Fold Change Shrinkage (apeglm) to clean up low-count noise
res_shrunk <- lfcShrink(dds, coef = "disease_JMML_vs_Normal", type = "apeglm")

# Convert results to a clean data frame
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("Symbol") %>%
  rename(log2FoldChange = log2FoldChange,
         pvalue         = pvalue,
         qvalue         = padj) %>% # DESeq2 uses padj (adjusted p-value)
  mutate(Status = case_when(
    qvalue < QVAL_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Upregulated",
    qvalue < QVAL_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Downregulated",
    TRUE                                                ~ "Not Significant"
  ))

print(table(res_df$Status, useNA = "ifany"))

# Save expression tables
write.csv(res_df, file.path(OUTPUT_DIR, "DEGs_all.csv"), row.names = FALSE)
write.csv(res_df %>% filter(Status != "Not Significant"),
          file.path(OUTPUT_DIR, "DEGs_significant.csv"), row.names = FALSE)


# ---- 5. Volcano Plot Generation ----
n_up  <- sum(res_df$Status == "Upregulated", na.rm = TRUE)
n_down <- sum(res_df$Status == "Downregulated", na.rm = TRUE)

genes_to_label <- c("TNF", "IL1B", "NFKB1", "NFKBIA", "PTPN11", "CSF2RA",
                    "JAK2", "STAT5A", "MARCKSL1", "ROR2")
top_hits <- res_df %>%
  filter(Status != "Not Significant") %>%
  arrange(qvalue) %>%
  head(10) %>%
  pull(Symbol)

label_df <- res_df %>%
  filter(Symbol %in% unique(c(genes_to_label, top_hits)),
         Status != "Not Significant")

p_volcano <- ggplot(res_df, aes(log2FoldChange, -log10(qvalue), color = Status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Upregulated"     = "firebrick",
                                "Downregulated"   = "steelblue",
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed") +
  geom_hline(yintercept = -log10(QVAL_CUTOFF), linetype = "dashed") +
  geom_text_repel(data = label_df, aes(label = Symbol),
                  color = "black", size = 3, max.overlaps = 20) +
  labs(title = paste0("DESeq2: JMML vs Normal (Up: ", n_up, ", Down: ", n_down, ")"),
       x = "log2 Fold Change",
       y = "-log10 qvalue") +
  theme_bw()

ggsave(file.path(OUTPUT_DIR, "volcano.png"), p_volcano, width = 8, height = 6, dpi = 200)


# ---- 6. Map to Entrez IDs ----
gene_map <- bitr(res_df$Symbol, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db, drop = TRUE)
res_entrez <- inner_join(res_df, gene_map, by = c("Symbol" = "SYMBOL"))

# Drop NAs from statistical filters
res_entrez <- res_entrez %>% filter(!is.na(qvalue))

# Significant genes for Over-Representation Analysis
sig_entrez <- res_entrez %>%
  filter(Status != "Not Significant") %>%
  pull(ENTREZID)

# Ranked list for GSEA
ranked <- res_entrez$log2FoldChange
names(ranked) <- res_entrez$ENTREZID
ranked <- sort(ranked[!duplicated(names(ranked))], decreasing = TRUE)


# ---- 7. KEGG Pathway Analysis ----
message("Running KEGG enrichment...")
kegg_ora <- enrichKEGG(gene         = sig_entrez,
                       organism     = "hsa",
                       pvalueCutoff = 0.05)

if (!is.null(kegg_ora) && nrow(as.data.frame(kegg_ora)) > 0) {
  kegg_ora <- setReadable(kegg_ora, org.Hs.eg.db, keyType = "ENTREZID")
  
  p_kegg <- dotplot(kegg_ora, showCategory = 15, title = "KEGG Pathways (DESeq2)") +
    theme_bw()
  ggsave(file.path(OUTPUT_DIR, "KEGG_dotplot.png"), p_kegg, width = 9, height = 7, dpi = 200)
  write.csv(as.data.frame(kegg_ora), file.path(OUTPUT_DIR, "KEGG_results.csv"), row.names = FALSE)
}


# ---- 8. Hallmark GSEA ----
message("Running Hallmark GSEA...")
hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, entrez_gene)

hallmark_gsea <- GSEA(geneList      = ranked,
                      TERM2GENE     = hallmarks,
                      minGSSize     = 10,
                      maxGSSize     = 500,
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      verbose       = FALSE)

if (!is.null(hallmark_gsea) && nrow(as.data.frame(hallmark_gsea)) > 0) {
  hallmark_gsea <- setReadable(hallmark_gsea, org.Hs.eg.db, keyType = "ENTREZID")
  hall_df <- as.data.frame(hallmark_gsea)
  
  top_hall <- hall_df %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(ID = gsub("HALLMARK_", "", ID),
           Direction = ifelse(NES > 0, "Up", "Down"),
           ID = factor(ID, levels = rev(ID)))
  
  p_hall <- ggplot(top_hall, aes(NES, ID, fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue")) +
    geom_vline(xintercept = 0) +
    labs(title = "Top Hallmark Pathways (DESeq2 GSEA)", x = "NES", y = NULL) +
    theme_bw()
  
  ggsave(file.path(OUTPUT_DIR, "Hallmark_barplot.png"), p_hall, width = 9, height = 7, dpi = 200)
  write.csv(hall_df, file.path(OUTPUT_DIR, "Hallmark_results.csv"), row.names = FALSE)
}


# ---- 9. Reactome Pathway Analysis ----
message("Running Reactome Enrichment Analysis...")
library(ReactomePA)

# Over-Representation Analysis (ORA) on significant DEGs
reactome_ora <- enrichPathway(gene         = sig_entrez,
                              organism     = "human",
                              pvalueCutoff = 0.05,
                              readable     = TRUE) # Automatically converts Entrez to Symbols

if (!is.null(reactome_ora) && nrow(as.data.frame(reactome_ora)) > 0) {
  p_reactome <- dotplot(reactome_ora, showCategory = 15, title = "Reactome Pathways (DESeq2)") +
    theme_bw()
  ggsave(file.path(OUTPUT_DIR, "Reactome_dotplot.png"), p_reactome, width = 10, height = 7, dpi = 200)
  write.csv(as.data.frame(reactome_ora), file.path(OUTPUT_DIR, "Reactome_results.csv"), row.names = FALSE)
}

# Gene Set Enrichment Analysis (GSEA) on all ranked genes
reactome_gsea <- gsePathway(geneList     = ranked,
                            organism     = "human",
                            minGSSize    = 10,
                            maxGSSize    = 500,
                            pvalueCutoff = 0.05,
                            verbose      = FALSE)

if (!is.null(reactome_gsea) && nrow(as.data.frame(reactome_gsea)) > 0) {
  reactome_gsea <- setReadable(reactome_gsea, org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(reactome_gsea), file.path(OUTPUT_DIR, "Reactome_GSEA_results.csv"), row.names = FALSE)
}



# ---- 10. Disease Ontology GSEA (Alternative to ORA) ----
message("Running Disease Ontology GSEA...")

dose_gsea <- gseDO(geneList     = ranked,
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose      = FALSE)

if (!is.null(dose_gsea) && nrow(as.data.frame(dose_gsea)) > 0) {
  # Convert Entrez IDs to readable Gene Symbols
  dose_gsea <- setReadable(dose_gsea, org.Hs.eg.db, keyType = "ENTREZID")
  dose_gsea_df <- as.data.frame(dose_gsea)
  
  # Clean up terms for a barplot or dotplot
  top_dose_gsea <- dose_gsea_df %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(Direction = ifelse(NES > 0, "Up", "Down"),
           Description = factor(Description, levels = rev(Description)))
  
  # Plot the Normalized Enrichment Scores (NES)
  p_dose_gsea <- ggplot(top_dose_gsea, aes(NES, Description, fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Top Disease Ontology Terms (DESeq2 GSEA)", x = "NES", y = NULL) +
    theme_bw()
  
  # Save the outputs
  ggsave(file.path(OUTPUT_DIR, "DOSE_GSEA_barplot.png"), p_dose_gsea, width = 10, height = 7, dpi = 200)
  write.csv(dose_gsea_df, file.path(OUTPUT_DIR, "DOSE_GSEA_results.csv"), row.names = FALSE)
  
  message("DOSE GSEA completed! Plot and table saved successfully.")
} else {
  message("No significant Disease Ontology terms found via GSEA either.")
}
