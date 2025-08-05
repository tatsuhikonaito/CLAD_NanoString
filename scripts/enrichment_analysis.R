# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(forcats)

if (requireNamespace("rstudioapi", quietly = TRUE)) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_dir <- dirname(script_path)
  print(script_dir)
}
base_dir <- dirname(script_dir)
setwd(base_dir)


OUT_dir <- file.path(base_dir, "geomxtools")
DEG_dir <- file.path(base_dir, "DEG")
enrichment_dir <- file.path(base_dir, "enrichment")

regions <- c("Epithelium", "Endothelium", "Bulk")
classes <- c("CLAD", "Non-CLAD", "Normal")
class_pairs = list(c("CLAD", "Non-CLAD"), c("CLAD", "Normal"), c("Non-CLAD", "Normal"))
class_pair_labels <- sapply(class_pairs, function(p) paste(p, collapse = " vs "))



region_colors <- setNames(RColorBrewer::brewer.pal(length(regions), "Set2"), regions)

# These were got mis_converted <- setdiff(sig_genes$Gene, gene_conversion$SYMBOL)
gene_conversion_df <- read.table("data/gene_entrez_conversion.txt", 
                                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_conversion_list <- as.list(setNames(gene_conversion_df$Entrez, gene_conversion_df$Gene))

mis_converted = c()


# Main loop
for (class_pair in class_pairs) {
  for (region in regions) {
    file_path <- file.path(DEG_dir, "results", paste("limma.results", paste(class_pair, collapse = "_"), region, "txt", sep = "."))
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    results_ <- read.table(file_path, sep = "\t", header = TRUE)
    
    gene_sets <- list(
      DEG_genes = results_ %>% filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>% pull(Gene),
      up_genes  = results_ %>% filter(logFC > 1 & adj.P.Val < 0.05) %>% pull(Gene),
      down_genes = results_ %>% filter(logFC < -1 & adj.P.Val < 0.05) %>% pull(Gene)
    )
    
    for (set_name in names(gene_sets)) {
      sig_genes <- gene_sets[[set_name]]
      sig_genes <- recode(sig_genes, !!!gene_conversion_list)
      
      sig_genes_entrez <- suppressMessages(
        bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      )
      mis_converted <- c(mis_converted, setdiff(sig_genes, sig_genes_entrez$SYMBOL))
      if (nrow(sig_genes_entrez) == 0) next
      
      # --- KEGG ---
      kegg_enrich <- enrichKEGG(
        gene = sig_genes_entrez$ENTREZID,
        organism = "hsa",
        pvalueCutoff = 1
      )
      
      if (!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
        enrich_df <- as.data.frame(kegg_enrich) %>%
          mutate(log10p_adj = -log10(p.adjust)) %>%
          arrange(desc(log10p_adj))
        
        file_prefix <- paste("limma.KEGG_enrichment", paste(class_pair, collapse = "_"), region, set_name, sep = ".")
        write.table(enrich_df, file.path(enrichment_dir, "results", paste0(file_prefix, ".txt")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
        # Optional: visualize with emapplot
        kegg_enrich <- pairwise_termsim(kegg_enrich)
        kegg_plot <- emapplot(kegg_enrich, showCategory = 15)
        
        ggsave(file.path(enrichment_dir, "plots", paste0(file_prefix, ".emapplot.pdf")),
               plot = kegg_plot, width = 8, height = 6)
      }
      
      # --- GO: Biological Process ---
      go_enrich <- enrichGO(
        gene = sig_genes_entrez$ENTREZID,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pvalueCutoff = 1
      )
      
      if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
        go_simplified <- simplify(go_enrich, cutoff = 0.7, by = "p.adjust", select_fun = min)
        
        file_prefix <- paste("limma.GO_BP_enrichment", paste(class_pair, collapse = "_"), region, set_name, sep = ".")
        
        write.table(as.data.frame(go_enrich), file.path(enrichment_dir, "results", paste0(file_prefix, ".txt")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        write.table(as.data.frame(go_simplified), file.path(enrichment_dir, "results", paste0(file_prefix, ".simplified.txt")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
        go_top <- as.data.frame(go_simplified) %>% filter(p.adjust < 0.05) %>% head(10)
        
        if (nrow(go_top) > 0) {
          go_plot <- ggplot(go_top, aes(x = fct_reorder(Description, log10(p.adjust)), y = -log10(p.adjust))) +
            geom_bar(stat = "identity", fill = "darkgreen") +
            coord_flip() +
            labs(x = "GO Biological Process", y = "−log10(adjusted P-value)",
                 title = paste0(class_pair[1], " vs ", class_pair[2], " — ", region, " — ", set_name)) +
            theme_bw(base_size = 10)
          
          ggsave(file.path(enrichment_dir, "plots", paste0(file_prefix, ".bar_plot.png")),
                 plot = go_plot, width = 8, height = 6, dpi = 300)
        }
      }
    }
  }
}
