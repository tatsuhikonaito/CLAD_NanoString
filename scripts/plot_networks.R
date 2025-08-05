# ==================== Load Libraries ====================
library(GeomxTools)
library(ComplexHeatmap)
library(RColorBrewer)
library(tibble)
library(circlize)
library(org.Hs.eg.db)
library(igraph)
library(ggraph)
library(tidygraph)
library(KEGGREST)
library(tidyr)
library(dplyr)
library(forcats)

# ==================== Paths and Setup ====================

if (requireNamespace("rstudioapi", quietly = TRUE)) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_dir <- dirname(script_path)
  print(script_dir)
}
base_dir <- dirname(script_dir)
setwd(base_dir)

plot_dir <- file.path(base_dir, "network", "plots")
OUT_dir <- file.path(base_dir, "geomxtools")
DEG_dir <- file.path(base_dir, "DEG")
enrichment_dir <- file.path(base_dir, "enrichment")
results_dir <- file.path(enrichment_dir, "results")

top_n_genes <- 20
regions <- c("Epithelium", "Endothelium", "Bulk")
classes <- c("CLAD", "Non-CLAD", "Normal")
class_pairs <- list(c("CLAD", "Non-CLAD"), c("CLAD", "Normal"), c("Non-CLAD", "Normal"))

class_colors <- c("CLAD" = "#C77CFF", "Non-CLAD" = "#00A0E0", "Normal" = "#E6C229")
region_colors <- setNames(brewer.pal(length(regions), "Set2"), regions)

# ==================== KEGG Genes for TNF / IL-17 ====================
tnf_kegg <- keggGet("hsa04668")[[1]]
tnf_entrez <- gsub("hsa:", "", tnf_kegg$GENE[seq(1, length(tnf_kegg$GENE), 2)])
tnf_genes <- mapIds(org.Hs.eg.db, keys = tnf_entrez, keytype = "ENTREZID", column = "SYMBOL")

il17_kegg <- keggGet("hsa04657")[[1]]
il17_entrez <- gsub("hsa:", "", il17_kegg$GENE[seq(1, length(il17_kegg$GENE), 2)])
il17_genes <- mapIds(org.Hs.eg.db, keys = il17_entrez, keytype = "ENTREZID", column = "SYMBOL")

pathway_genes <- union(tnf_genes, il17_genes)

# ==================== Load Expression ====================
target_demoData <- readRDS(file.path(OUT_dir, "target_demoData.norm.rds"))
exp_genes <- assayDataElement(target_demoData, "tmm_cpm")
metadata <- as.data.frame(pData(target_demoData)) %>% filter(region %in% regions)
exp_genes <- exp_genes[, colnames(exp_genes) %in% rownames(metadata)]
metadata <- metadata[colnames(exp_genes), ]
metadata$group <- paste(metadata$region, metadata$class, sep = "_")

# Normalize within regions (z-score)
exp_genes_norm <- matrix(NA, nrow = nrow(exp_genes), ncol = ncol(exp_genes),
                         dimnames = list(rownames(exp_genes), colnames(exp_genes)))
for (r in regions) {
  region_samples <- rownames(metadata)[metadata$region == r]
  mat <- exp_genes[, region_samples, drop = FALSE]
  exp_genes_norm[, region_samples] <- t(scale(t(mat)))
}

df_long <- as.data.frame(t(exp_genes_norm)) %>%
  mutate(individual = metadata$individual,
         group = metadata$group)

df_indiv_avg <- df_long %>%
  group_by(individual, group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

grouped_expression <- df_indiv_avg %>%
  unite("sample_id", individual, group, sep = "__") %>%
  column_to_rownames("sample_id") %>%
  t()

sample_info <- data.frame(
  sample_id = colnames(grouped_expression),
  individual = sapply(strsplit(colnames(grouped_expression), "__"), `[`, 1),
  group = sapply(strsplit(colnames(grouped_expression), "__"), `[`, 2),
  stringsAsFactors = FALSE
)
sample_info$region <- sapply(strsplit(sample_info$group, "_"), `[`, 1)
sample_info$class <- sapply(strsplit(sample_info$group, "_"), `[`, 2)
sample_info$class <- factor(sample_info$class, levels = c("CLAD", "Non-CLAD", "Normal"))


class_pair = c("CLAD", "Non-CLAD")
class_tag <- paste(class_pair, collapse = "_")
set_name = "up_genes"
region1 <- "Epithelium"
region2 <- "Endothelium"

# Function to load and filter one region
load_and_filter_enrich <- function(region) {
  file_prefix <- paste("limma.KEGG_enrichment", class_tag, region, set_name, sep = ".")
  file_path <- file.path(results_dir, paste0(file_prefix, ".txt"))
  
  read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    filter(subcategory %in% c("Immune system", "Signal transduction")) %>%
    filter(p.adjust < 0.05) %>%
    mutate(log10p = -log10(p.adjust)) %>%
    arrange(p.adjust) %>%
    mutate(Description = fct_reorder(Description, log10p))
}

# Load and filter each
filtered1 <- load_and_filter_enrich(region1)
filtered2 <- load_and_filter_enrich(region2)

# Get shared KEGG terms
shared_descriptions <- intersect(filtered1$Description, filtered2$Description)
# This gets "IL-17 signaling pathway" "TNF signaling pathway"
# Accordingly, we focuses on these pathways in the following.

for (region in regions) {
  pair_name <- paste(class_pair, collapse = "_")
  file_path <- file.path(DEG_dir, "results", paste("limma.results", pair_name, region, "txt", sep = "."))
  
  deg_df <- read.table(file_path, sep = "\t", header = TRUE)
  deg_df$gene <- deg_df$Gene
  region_samples <- sample_info$sample_id[sample_info$region == region]
  
  # -------------------- [3] Pathway Network Plot --------------------
  deg_pathway <- deg_df %>% filter(gene %in% pathway_genes)
  lfc_thresh <- 0.5
  deg_pathway <- deg_pathway %>% filter(abs(logFC) >= lfc_thresh)
  
  genes_for_net <- deg_pathway$gene[deg_pathway$gene %in% rownames(grouped_expression)]
  if (length(genes_for_net) >= 3) {
    expr_sub <- grouped_expression[genes_for_net, region_samples, drop = FALSE]
    expr_sub <- expr_sub[rowSums(is.na(expr_sub)) == 0, , drop = FALSE]
    
    if (nrow(expr_sub) >= 3) {
      cor_mat <- cor(t(expr_sub), method = "pearson", use = "pairwise.complete.obs")
      cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
      
      cor_thresh <- 0.7
      edge_idx <- which(abs(cor_mat) > cor_thresh, arr.ind = TRUE)
      
      if (nrow(edge_idx) > 1000) {
        message("Too many edges: skipping network plot.")
        next
      }
      
      # Build edges with sign info
      edges <- data.frame(
        from = rownames(cor_mat)[edge_idx[, 1]],
        to = colnames(cor_mat)[edge_idx[, 2]],
        weight = abs(cor_mat[edge_idx]),
        sign = ifelse(cor_mat[edge_idx] > 0, "positive", "negative"),
        linetype = ifelse(cor_mat[edge_idx] > 0, "solid", "dashed")
      ) %>%
        filter(!is.na(from) & !is.na(to) & from != to) %>%
        dplyr::select(from, to, weight, sign, linetype)
      
      nodes <- deg_pathway %>%
        filter(gene %in% union(edges$from, edges$to)) %>%
        mutate(logFC_abs = abs(logFC),
               Pathway = case_when(
                 gene %in% tnf_genes & gene %in% il17_genes ~ "Both",
                 gene %in% tnf_genes ~ "TNF",
                 gene %in% il17_genes ~ "IL-17",
                 TRUE ~ "Other"
               )) %>%
        distinct(gene, .keep_all = TRUE) %>%
        rename(name = gene)
      
      nodes$PathwayList <- lapply(nodes$name, function(gene) {
        pathways <- character(0)
        if (gene %in% il17_genes) pathways <- c(pathways, "IL-17")
        if (gene %in% tnf_genes) pathways <- c(pathways, "TNF")
        return(pathways)
      })
      
      edges$layout_weight <- edges$weight
      g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
      E(g)$weight <- edges$layout_weight
      g_tbl <- as_tbl_graph(g)
      
      png(file.path(plot_dir, paste0("network.", region, ".", pair_name, ".png")), width = 1200, height = 960, res = 250)
      print(
        ggraph(g_tbl, layout = "fr") +
          geom_edge_link(
            aes(edge_alpha = weight, linetype = linetype),
            edge_width = 0.6,
            color = "grey30",
            show.legend = FALSE
          ) +
          geom_node_point(aes(size = logFC_abs, fill = Pathway), 
                          color = "black", shape = 21, stroke = 0.3) +
          geom_node_text(aes(label = name), repel = TRUE, size = 2.5, segment.size = 0.2) +
          
          # LogFC-based node size legend
          scale_size_continuous(
            name = expression("|" * log[2] * "FC" * "|"),
            range = c(1.5, 7),
            breaks = c(1, 2, 3)
          ) +
          
          # Pathway-based fill color legend
          scale_fill_manual(values = c(
            "TNF"   = "#EF6C42",
            "IL-17" = "#B27ACD",
            "Both"  = "#D64976"
          )) +
          
          guides(
            fill = guide_legend(
              override.aes = list(size = 4),
              title = "Pathway",
              order = 1  # ensures Pathway comes first
            ),
            size = guide_legend(
              override.aes = list(shape = 21, fill = "grey50", color = "black"),
              title = expression("|" * log[2] * "FC" * "|"),
              order = 2  # ensures log2FC comes after
            )
          )
        +
          
          scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
          theme_void()
      )
      
      dev.off()
    }
  }
}

