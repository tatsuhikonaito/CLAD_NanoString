
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)


base_dir <- "/Users/tatsuhiko/Documents/Projects/lung_transplantation"
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


for (class_pair in class_pairs) {
  for (region in regions) {
    file_path <- file.path(DEG_dir, "results", paste("limma.results", paste(class_pair, collapse = "_"), region, "txt", sep = "."))
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    results_ <- read.table(file_path, sep = "\t", header = TRUE)
    
    # Create gene sets
    gene_sets <- list(
      DEG_genes = results_ %>% filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>% pull(Gene),
      up_genes  = results_ %>% filter(logFC > 1 & adj.P.Val < 0.05) %>% pull(Gene),
      down_genes = results_ %>% filter(logFC < -1 & adj.P.Val < 0.05) %>% pull(Gene)
    )
    
    for (set_name in names(gene_sets)) {
      sig_genes <- gene_sets[[set_name]]
      
      # Recode gene symbols (if any corrections apply)
      sig_genes <- recode(sig_genes, !!!gene_conversion_list)
      
      # Convert to Entrez IDs
      sig_genes_entrez <- suppressMessages(
        bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      )
      
      mis_converted <- c(mis_converted, setdiff(sig_genes, sig_genes_entrez$SYMBOL))
      
      # Skip if conversion failed
      if (nrow(sig_genes_entrez) == 0) next
      
      ### --- KEGG ---
      kegg_enrich <- enrichKEGG(
        gene = sig_genes_entrez$ENTREZID,
        organism = "hsa",
        pvalueCutoff = 1  # <- include all results
      )
      
      if (!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
        enrich_df <- as.data.frame(kegg_enrich) %>%
          mutate(log10p_adj = -log10(p.adjust)) %>%
          arrange(desc(log10p_adj))
        
        file_prefix <- paste("limma.KEGG_enrichment", paste(class_pair, collapse = "_"), region, set_name, sep = ".")
        
        # ✅ Save all enrichment results
        write.table(enrich_df, file.path(enrichment_dir, "results", paste0(file_prefix, ".txt")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
        # ✅ Filter for plotting (p.adjust < 0.05)
        enrich_top <- enrich_df %>% filter(p.adjust < 0.05) %>% head(10)
        
        if (nrow(enrich_top) > 0) {
          kegg_plot <- ggplot(enrich_top, aes(x = reorder(Description, log10p_adj), y = log10p_adj)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() +
            labs(x = "KEGG pathway", y = "−log10(adjusted P-value)",
                 title = paste0(class_pair[1], " vs ", class_pair[2], " — ", region, " — ", set_name)) +
            theme_bw(base_size = 10)
          
          ggsave(file.path(enrichment_dir, "plots", paste0(file_prefix, ".bar_plot.png")),
                 plot = kegg_plot, width = 8, height = 6, dpi = 300)
        }
      }
      
      ### --- GO BP ---
      go_enrich <- enrichGO(
        gene = sig_genes_entrez$ENTREZID,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pvalueCutoff = 1  # <- changed from 0.05 to 1 to keep all results
      )
      
      if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
        go_df <- as.data.frame(go_enrich) %>%
          mutate(log10p_adj = -log10(p.adjust)) %>%
          arrange(desc(log10p_adj))
        
        file_prefix <- paste("limma.GO_BP_enrichment", paste(class_pair, collapse = "_"), region, set_name, sep = ".")
        
        # ✅ Save all results
        write.table(go_df, file.path(enrichment_dir, "results", paste0(file_prefix, ".txt")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
        # ✅ Filter for plotting only
        go_top <- go_df %>% filter(p.adjust < 0.05) %>% head(10)
        
        if (nrow(go_top) > 0) {
          go_plot <- ggplot(go_top, aes(x = reorder(Description, log10p_adj), y = log10p_adj)) +
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

class_pair = c("CLAD", "Non-CLAD")

gene_sets <- c("up_genes", "down_genes")


for (region in regions) {
  region_results <- list()
  sig_terms <- c()
  
  # Step 1: Load results and collect significant GO Descriptions
  for (set_name in gene_sets) {
    file_prefix <- paste("limma.KEGG_enrichment", paste(class_pair, collapse = "_"), region, set_name, sep = ".")
    result_path <- file.path(enrichment_dir, "results", paste0(file_prefix, ".txt"))
    
    if (file.exists(result_path)) {
      df <- read.table(result_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
        mutate(set_name = set_name)
      
      region_results[[set_name]] <- df
      
      # Collect significant GO terms
      sig_terms <- union(sig_terms, df %>% filter(p.adjust < 0.05) %>% pull(Description))
    }
  }
  
  if (length(sig_terms) == 0) next
  
  combined_df <- bind_rows(region_results) %>%
    filter(Description %in% sig_terms) %>%
    dplyr::select(Description, set_name, p.adjust) %>%
    mutate(
      log2p = -log2(p.adjust),
      set_name = factor(set_name, levels = c("down_genes", "up_genes"))  # down below, up above
    )
  
  plot <- ggplot(combined_df, aes(x = Description, y = set_name)) +
    geom_point(aes(size = log2p, fill = set_name), shape = 21, color = "black") +
    scale_size_continuous(name = expression(-log[2](adjusted~P~value))) +
    scale_fill_manual(values = c("up_genes" = "tomato", "down_genes" = "steelblue")) +
    labs(
      title = paste("KEGG pathway enrichment in", region),
      x = "Pathway",
      y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.position = "right"
    )
  
  ggsave(file.path(enrichment_dir, "plots", paste0("KEGG_dotplot_", paste(class_pair, collapse = "_"), ". ", region, ".png")),
         plot = plot, width = max(8, length(sig_terms) * 0.3), height = 4, dpi = 300)
}





class_pairs_target <- list(c("CLAD", "Normal"), c("Non-CLAD", "Normal"))

for (region in regions) {
  results_list <- list()
  sig_terms <- c()
  
  for (class_pair in class_pairs_target) {
    class_label <- paste(class_pair, collapse = "_")
    
    for (set_name in gene_sets) {
      file_prefix <- paste("limma.KEGG_enrichment", class_label, region, set_name, sep = ".")
      result_path <- file.path(enrichment_dir, "results", paste0(file_prefix, ".txt"))
      
      if (file.exists(result_path)) {
        df <- read.table(result_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
          mutate(
            set_name = set_name,
            class_label = class_label,
            combo = paste(set_name, class_label, sep = " | ")
          )
        
        results_list[[paste(region, class_label, set_name, sep = "_")]] <- df
        
        # Collect significant GO terms
        sig_terms <- union(sig_terms, df %>% filter(p.adjust < 0.05) %>% pull(Description))
      }
    }
  }
  
  if (length(sig_terms) == 0) next
  
  # Combine and filter to significant terms
  combined_df <- bind_rows(results_list) %>%
    filter(Description %in% sig_terms) %>%
    mutate(
      log2p = -log2(p.adjust),
      combo = factor(combo, levels = rev(unique(combo)))  # ensures up above, down below
    ) %>%
    dplyr::select(Description, combo, set_name, log2p)
  
  # Dot plot
  plot <- ggplot(combined_df, aes(x = Description, y = combo)) +
    geom_point(aes(size = log2p, fill = set_name), shape = 21, color = "black") +
    scale_size_continuous(name = expression(-log[2](adjusted~P~value))) +
    scale_fill_manual(values = c("up_genes" = "tomato", "down_genes" = "steelblue")) +
    labs(
      title = paste("KEGG pathway enrichment in", region),
      x = "Pathway",
      y = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  # Save the plot
  ggsave(file.path(enrichment_dir, "plots", paste0("KEGG_dotplot_classpairs_", region, ".png")),
         plot = plot, width = max(10, length(sig_terms) * 0.3), height = 5, dpi = 300)
}