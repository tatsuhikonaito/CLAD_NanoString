
library(dplyr)
library(ggplot2)

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
classes <- c("Non-CLAD", "CLAD", "Normal")
class_pairs <- list(c("CLAD", "Non-CLAD"), c("Non-CLAD", "Normal"))
class_pair_labels <- sapply(class_pairs, function(p) paste(p, collapse = " vs "))

# --- Count significant genes ---
sig_gene_counts <- data.frame()

for (region in regions) {
  for (class_pair in class_pairs) {
    file_path <- file.path(res_dir, paste("limma.results", paste(class_pair, collapse = "_"), region, "txt", sep = "."))
    if (file.exists(file_path)) {
      results_ <- read.table(file_path, sep = "\t", header = TRUE)
      sig_genes <- results_ %>%
        filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
        pull(Gene)
      
      sig_gene_counts <- rbind(sig_gene_counts, data.frame(
        region = region,
        class1 = class_pair[1],
        class2 = class_pair[2],
        n_sig_genes = length(sig_genes)
      ))
    } else {
      warning(paste("File not found:", file_path))
    }
  }
}

# --- Add x-axis label column ---
sig_gene_counts <- sig_gene_counts %>%
  mutate(
    class_pair = paste(class1, class2, sep = " vs "),
    x_label = paste(class_pair, region, sep = " | ")
  )

# --- Define full x-axis label levels with spacers ---
spacer_labels <- paste0("spacer_", seq_len(length(class_pair_labels) - 1))
x_levels_spaced <- c()

for (i in seq_along(class_pair_labels)) {
  pair <- class_pair_labels[i]
  group_levels <- paste(pair, regions, sep = " | ")
  x_levels_spaced <- c(x_levels_spaced, group_levels)
  if (i < length(class_pair_labels)) {
    x_levels_spaced <- c(x_levels_spaced, spacer_labels[i])
  }
}
x_levels_spaced <- c("padding_left", x_levels_spaced, "padding_right")

# --- Prepare plot data ---
sig_gene_counts_plot <- sig_gene_counts %>%
  mutate(x_label = factor(x_label, levels = x_levels_spaced))

spacing_rows <- data.frame(
  class1 = NA, class2 = NA, region = NA,
  n_sig_genes = NA, class_pair = NA,
  x_label = factor(spacer_labels, levels = x_levels_spaced)
)

padding_rows <- data.frame(
  class1 = NA, class2 = NA, region = NA,
  n_sig_genes = NA, class_pair = NA,
  x_label = factor(c("padding_left", "padding_right"), levels = x_levels_spaced)
)

plot_df <- bind_rows(sig_gene_counts_plot, spacing_rows, padding_rows)

# --- Axis labels for regions only (or blank for spacers) ---
x_axis_labels <- setNames(
  sapply(x_levels_spaced, function(x) {
    if (grepl("^spacer_|^padding_", x)) return("")
    strsplit(x, " \\| ")[[1]][2]
  }),
  x_levels_spaced
)

# --- Region factor order ---
plot_df$region <- factor(plot_df$region, levels = regions)

# --- Dynamically recompute class pair labels and center positions ---
used_class_pairs <- unique(na.omit(plot_df$class_pair))
center_positions <- sapply(used_class_pairs, function(pair) {
  group_labels <- paste(pair, regions, sep = " | ")
  idx <- which(levels(plot_df$x_label) %in% group_labels)
  mean(idx)
})

# --- Build the plot ---
plt <- ggplot(plot_df, aes(x = x_label, y = n_sig_genes, fill = region)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", na.rm = TRUE) +
  scale_x_discrete(labels = rep("", length(x_levels_spaced))) +
  scale_fill_manual(values = region_colors, na.translate = FALSE) +
  labs(x = NULL, y = "Differentially expressed genes") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "right",
    plot.margin = margin(20, 20, 40, 20)
  ) +
  annotate("text",
           x = center_positions,
           y = -max(plot_df$n_sig_genes, na.rm = TRUE) * 0.05,
           label = used_class_pairs,
           angle = 0, hjust = 0.5, size = 4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4)

# Save
ggsave(filename = file.path(plot_dir, "DEG_genes.bar_plot.png"),
       plot = plt, width = 7, height = 5, dpi = 300)

## If we show three pairs
#ggsave(filename = file.path(plot_dir, paste("DEG_genes.bar_plot.png", sep=".")), 
#       plot = plt, width = 9, height = 5, dpi = 300)


region_1 <- "Epithelium"
region_2 <- "Endothelium"


class_pair <- c("CLAD", "Non-CLAD")
  
results_1 <- read.table(file.path(res_dir, paste("limma.results", paste(class_pair, collapse = "_"), region_1, "txt", sep=".")),
                        sep = "\t", header = T)
sig_genes_1 <- results_1 %>% filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>% pull(Gene)


results_2 <- read.table(file.path(res_dir, paste("limma.results", paste(class_pair, collapse = "_"), region_2, "txt", sep=".")),
                        sep = "\t", header = T)
sig_genes_2 <- results_2 %>% filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>% pull(Gene)



library(VennDiagram)

# Example vectors (replace with your actual gene lists)
# sig_genes_1 <- c("IFNG", "STAT1", "IRF1")
# sig_genes_2 <- c("STAT1", "IRF9", "MX1")

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = list(
    "Set 1" = sig_genes_1,
    "Set 2" = sig_genes_2
  ),
  filename = NULL,  # Prevent automatic file saving
  fill = c("skyblue", "lightpink"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.05,
  margin = 0.1
)

# Draw the plot
grid.newpage()
grid.draw(venn.plot)



