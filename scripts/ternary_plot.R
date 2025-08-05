library(GeomxTools)
library(ggtern)
library(tidyverse)

if (requireNamespace("rstudioapi", quietly = TRUE)) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_dir <- dirname(script_path)
  print(script_dir)
}
base_dir <- dirname(script_dir)
setwd(base_dir)

plot_dir <- file.path(base_dir, "DEG", "plots")
OUT_dir <- file.path(base_dir, "geomxtools")

# Define regions and classes
regions <- c("Epithelium", "Endothelium", "Bulk")
classes <- c("CLAD", "Non-CLAD", "Normal")

# Load expression and metadata
target_demoData <- readRDS(file.path(OUT_dir, "target_demoData.norm.rds"))
expr <- assayDataElement(target_demoData, "tmm_cpm")  # already on linear scale
metadata <- as.data.frame(pData(target_demoData)) %>%
  filter(region %in% regions)
expr <- expr[, colnames(expr) %in% rownames(metadata)]
metadata <- metadata[colnames(expr), ]

# Define group: region + class
metadata$group <- paste(metadata$region, metadata$class, sep = "_")
group_levels <- unlist(lapply(regions, function(r) paste(r, classes, sep = "_")))
metadata$group <- factor(metadata$group, levels = group_levels)

# Average per individual
df_long <- as.data.frame(t(expr)) %>%
  mutate(individual = metadata$individual,
         group = metadata$group)
df_indiv_avg <- df_long %>%
  group_by(individual, group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

# Then average per group
grouped_expression <- df_indiv_avg %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("group") %>%
  t()

# Plot per region
for (region in regions) {
  group_names <- paste(region, classes, sep = "_")
  if (!all(group_names %in% colnames(grouped_expression))) next
  
  # Normalize expression across the 3 classes
  df_tern <- grouped_expression[, group_names, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    filter(rowSums(is.na(across(all_of(group_names)))) == 0) %>%
    mutate(total = rowSums(across(all_of(group_names)))) %>%
    filter(total > 0) %>%
    mutate(across(all_of(group_names), ~ .x / total)) %>%
    rename(CLAD = paste0(region, "_CLAD"),
           NonCLAD = paste0(region, "_Non-CLAD"),
           Normal = paste0(region, "_Normal"))
  
  # --- Load DEGs for all comparisons ---
  deg_files <- list(
    CLAD_NonCLAD = file.path(base_dir, "DEG", "results", paste("limma.results.CLAD_Non-CLAD", region, "txt", sep = ".")),
    CLAD_Normal = file.path(base_dir, "DEG", "results", paste("limma.results.CLAD_Normal", region, "txt", sep = ".")),
    NonCLAD_Normal = file.path(base_dir, "DEG", "results", paste("limma.results.Non-CLAD_Normal", region, "txt", sep = "."))
  )
  
  deg_all_combined <- lapply(names(deg_files), function(name) {
    read_tsv(deg_files[[name]], show_col_types = FALSE) %>%
      filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
      mutate(pair = name)
  }) %>%
    bind_rows()
  
  # --- Assign up_class regardless of direction ---
  deg_with_class <- deg_all_combined %>%
    mutate(
      up_class = case_when(
        pair == "CLAD_NonCLAD" & logFC > 0 ~ "CLAD",
        pair == "CLAD_NonCLAD" & logFC < 0 ~ "Non-CLAD",
        pair == "CLAD_Normal"   & logFC > 0 ~ "CLAD",
        pair == "CLAD_Normal"   & logFC < 0 ~ "Normal",
        pair == "NonCLAD_Normal" & logFC > 0 ~ "Non-CLAD",
        pair == "NonCLAD_Normal" & logFC < 0 ~ "Normal",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(up_class))  # remove ambiguous ones
  
  # Use first up_class if gene appears in multiple pairs
  deg_up_final <- deg_with_class %>%
    group_by(Gene) %>%
    slice_min(order_by = P.Value, n = 1) %>%
    ungroup()
  
  # --- Merge into df_tern and assign color ---
  df_tern <- df_tern %>%
    left_join(deg_up_final, by = "Gene") %>%
    rowwise() %>%
    mutate(
      max_class = names(which.max(c(CLAD = CLAD, `Non-CLAD` = NonCLAD, Normal = Normal))),
      color = case_when(
        up_class == max_class ~ case_when(
          up_class == "CLAD" ~ "#C77CFF",   # Soft Magenta
          up_class == "Non-CLAD" ~ "#00A0E0",       # Cool Cyan
          up_class == "Normal" ~ "#E6C229"      # Muted Yellow (Improved visibility)
        ),
        TRUE ~ "gray"
      )
    ) %>%
    ungroup()
  
  # --- Top 20 genes per comparison (irrespective of direction) ---
  deg_top20_all <- deg_all_combined %>%
    group_by(pair) %>%
    arrange(P.Value) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  df_top_labels <- df_tern %>%
    semi_join(deg_top20_all, by = "Gene")
  
  # --- Plot ---
  p <- ggtern(df_tern, aes(x = NonCLAD, y = CLAD, z = Normal)) +
    geom_point(aes(color = color), size = 0.3, alpha = 0.6) +
    geom_text(
      data = df_top_labels,
      aes(label = Gene),
      size = 2.5,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_color_identity() +
    coord_tern(Tlim = 1, Llim = 1, Rlim = 1) +
    labs(
      T = "CLAD", L = "Non-CLAD", R = "Normal"
    ) +
    theme_void() +
    theme(
      tern.axis.title.T = element_text(size = 12, face = "bold"),
      tern.axis.title.L = element_text(size = 12, face = "bold"),
      tern.axis.title.R = element_text(size = 12, face = "bold"),
      tern.axis.text.T  = element_blank(),
      tern.axis.text.L  = element_blank(),
      tern.axis.text.R  = element_blank(),
      tern.axis.ticks.length.major = unit(0, "pt"),
      tern.axis.ticks.length.minor = unit(0, "pt"),
      tern.axis.arrow.show = FALSE,
      plot.margin = margin(20, 20, 20, 20)  # <- key to prevent text clipping
    )
  
  
  ggsave(file.path(plot_dir, paste0("ternary_", region, ".colored.png")),
         plot = p, width = 6, height = 5, dpi = 200)
  
}
