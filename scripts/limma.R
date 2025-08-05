library(edgeR)
library(limma)
library(ggrepel)
library(GeomxTools)

if (requireNamespace("rstudioapi", quietly = TRUE)) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_dir <- dirname(script_path)
  print(script_dir)
}
base_dir <- dirname(script_dir)
setwd(base_dir)

OUT_dir <- file.path(base_dir, "geomxtools")
res_dir <- file.path(base_dir, "DEG", "results")
plot_dir <- file.path(base_dir, "DEG", "plots")

regions <- c("Epithelium", "Endothelium", "Bulk")
classes <- c("CLAD", "Non-CLAD", "Normal")
class_pairs = list(c("CLAD", "Non-CLAD"), c("CLAD", "Normal"), c("Non-CLAD", "Normal"))


target_demoData <- readRDS(file.path(OUT_dir, "target_demoData.norm.rds"))

class_pair = class_pairs[[1]]
region = regions[1]

include_unique <- "unique"
include_unique <- "all"


for (class_pair in class_pairs){
  
  pData(target_demoData)$testClass <-
    factor(pData(target_demoData)$class, class_pair)
  
  for(region in regions) {
    include <- (pData(target_demoData)$region == region) & (pData(target_demoData)$testClass %in% class_pair)
    
    counts_matrix <- assayDataElement(object = target_demoData, elt = "exprs")[, include]
    counts_matrix <- counts_matrix[rownames(counts_matrix) != "NegProbe-WTX", ]
    
    # Extract sample metadata from the colData slotn
    metadata <- as.data.frame(pData(target_demoData))[include, ]
    metadata$testClass <- relevel(metadata$testClass, ref = class_pair[2])

    # Create a DGEList object and normalize
    dge <- DGEList(counts = counts_matrix)

    dge <- calcNormFactors(dge, method = "TMM")
    
    # Step 2: Create the design matrix
    design <- model.matrix(~ testClass, data = metadata)
    
    # Step 3: First voom transformation (without correlation)
    v1 <- voom(dge, design, plot = TRUE)
    
    # Step 4: Estimate correlation from voom-transformed data
    corfit <- duplicateCorrelation(v1, design, block = metadata$individual)
    cat("Estimated consensus correlation (1st pass):", corfit$consensus.correlation, "\n")
    
    # Step 5: Re-run voom, now incorporating correlation and blocking
    v2 <- voom(dge, design, plot = TRUE,
               block = metadata$individual,
               correlation = corfit$consensus.correlation)
    
    # Step 6: Re-estimate correlation (optional but recommended for precision)
    corfit2 <- duplicateCorrelation(v2, design, block = metadata$individual)
    cat("Updated consensus correlation (2nd pass):", corfit2$consensus.correlation, "\n")
    
    # Step 7: Fit the linear model using updated correlation
    fit <- lmFit(v2, design,
                 block = metadata$individual,
                 correlation = corfit2$consensus.correlation)
    fit <- eBayes(fit)
    
    # Extract the results; adjust the coefficient index or name as needed
    results <- topTable(fit, coef = paste0("testClass", class_pair[1]), number = Inf)
    results <- cbind(Gene = rownames(results), results)
    rownames(results) <- NULL
    head(results)

    
    write.table(results, file.path(res_dir, paste("limma.results", paste(class_pair, collapse = "_"), region, "txt", sep=".")),
                quote = F, sep = "\t", col.names = T, row.names = F)
    
    
    # Categorize Results based on P-value & FDR for plotting
    results$Color <- "NS or FC < 1"
    results$Color[results$adj.P.Val  < 0.05] <- "Adj.P < 0.05"
    results$Color[results$adj.P.Val < 0.001] <- "Adj.P < 0.001"
    results$Color[abs(results$logFC) < 1] <- "NS or FC < 1"
    results$Color <- factor(results$Color,
                            levels = c("NS or FC < 1", 
                                       "Adj.P < 0.05", "Adj.P < 0.001"))
    
    # pick top genes for either side of volcano to label
    # order genes for convenience:
    results$invert_P <- (-log10(results$P.Value)) * sign(results$logFC)
    
    top_g <- c()
    for(cond in regions) {
      top_g <- c(top_g,
                 results[, 'Gene'][
                   order(results[, 'invert_P'], decreasing = TRUE)[1:15]],
                 results[, 'Gene'][
                   order(results[, 'invert_P'], decreasing = FALSE)[1:15]])
    }
    top_g <- unique(top_g) 
    
    # Graph results
    # Compute maximum absolute log2 fold change to set symmetric limits on the x-axis
    max_abs <- max(abs(results$logFC), na.rm = TRUE) * 1.05
    
    plt <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val),
                               color = Color, label = Gene)) +
      geom_vline(xintercept = 0, color = "black", linetype = "solid") +  # center line at 0
      geom_vline(xintercept = c(1, -1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_point() +
      labs(x = bquote("Upregulated in" ~ .(class_pair[2]) ~ "<-" ~ log[2](FC) ~ "-> Upregulated in" ~ .(class_pair[1])),
           y = expression(-log[10]("adj." * italic(P))),
           title = paste0(class_pair[1], " vs ", class_pair[2], ", ", region),
           color = "Significance") +
      scale_color_manual(values = c(`Adj.P < 0.001` = "dodgerblue",
                                    `Adj.P < 0.05` = "lightblue",
                                    `NS or FC < 1` = "gray"),
                         guide = guide_legend(override.aes = list(size = 4))) +
      scale_x_continuous(limits = c(-max_abs, max_abs), expand = expansion(mult = c(0, 0.05))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      geom_text_repel(data = subset(results, Gene %in% top_g & adj.P.Val < 0.001),
                      size = 4, point.padding = 0.15, color = "black",
                      min.segment.length = 0.1, box.padding = 0.2, lwd = 2,
                      max.overlaps = 50) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))
    ggsave(filename = file.path(plot_dir, paste("limma.results", paste(class_pair, collapse = "_"), region, "volcano_plot.png", sep=".")), 
           plot = plt, width = 6, height = 6, dpi = 300)
    
  }
}

genes <- scan("data/IFN_genes.txt", "")


for (class_pair in class_pairs){
  
  pData(target_demoData)$testClass <-
    factor(pData(target_demoData)$class, class_pair)
  
  for(region in regions) {

    
    results <- read.table(file.path(res_dir, paste("limma.results", paste(class_pair, collapse = "_"), region, "txt", sep=".")),
                          sep = "\t", header = T)
    results <- results %>% filter(Gene %in% genes)
    
    # Categorize Results based on P-value & FDR for plotting
    results$Color <- "NS or FC < 1"
    #    results$Color[results$pvalue < 0.05] <- "P < 0.05"
    results$Color[results$adj.P.Val  < 0.05] <- "Adj.P < 0.05"
    results$Color[results$adj.P.Val < 0.001] <- "Adj.P < 0.001"
    results$Color[abs(results$logFC) < 1] <- "NS or FC < 1"
    results$Color <- factor(results$Color,
                            levels = c("NS or FC < 1", 
                                       "Adj.P < 0.05", "Adj.P < 0.001"))
    
    # pick top genes for either side of volcano to label
    # order genes for convenience:
    results$invert_P <- (-log10(results$P.Value)) * sign(results$logFC)
    
    top_g <- c()
    for(cond in regions) {
      top_g <- c(top_g,
                 results[, 'Gene'][
                   order(results[, 'invert_P'], decreasing = TRUE)[1:15]],
                 results[, 'Gene'][
                   order(results[, 'invert_P'], decreasing = FALSE)[1:15]])
    }
    top_g <- unique(top_g) 
    top_g <- genes
    
    # Graph results
    # Compute maximum absolute log2 fold change to set symmetric limits on the x-axis
    max_abs <- max(abs(results$logFC), na.rm = TRUE) * 1.05
    
    plt <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val),
                               color = Color, label = Gene)) +
      geom_vline(xintercept = 0, color = "black", linetype = "solid") +  # center line at 0
      geom_vline(xintercept = c(1, -1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_point() +
      labs(x = bquote("Upregulated in" ~ .(class_pair[2]) ~ "<-" ~ log[2](FC) ~ "-> Upregulated in" ~ .(class_pair[1])),
           y = expression(-log[10]("adj." * italic(P))),
           title = paste0(class_pair[1], " vs ", class_pair[2], ", ", region),
           color = "Significance") +
      scale_color_manual(values = c(`Adj.P < 0.001` = "dodgerblue",
                                    `Adj.P < 0.05` = "lightblue",
                                    `NS or FC < 1` = "gray"),
                         guide = guide_legend(override.aes = list(size = 4))) +
      scale_x_continuous(limits = c(-max_abs, max_abs), expand = expansion(mult = c(0, 0.05))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
#      geom_text_repel(data = subset(results, Gene %in% top_g & adj.P.Val < 0.001),
      geom_text_repel(data = subset(results, Gene %in% top_g),
                      size = 4, point.padding = 0.15, color = "black",
                      min.segment.length = 0.1, box.padding = 0.2, lwd = 2,
                      max.overlaps = 50) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))
    ggsave(filename = file.path(plot_dir, paste("limma.results", paste(class_pair, collapse = "_"), region, "IFN_genes.volcano_plot.png", sep=".")), 
           plot = plt, width = 6, height = 6, dpi = 300)
    
  }
}
