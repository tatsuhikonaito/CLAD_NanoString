# This script is based on the official GeoMx RNA-NGS analysis workflow:
# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(htmlwidgets)
library(webshot)

if (requireNamespace("rstudioapi", quietly = TRUE)) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_dir <- dirname(script_path)
  print(script_dir)
}
base_dir <- dirname(script_dir)
setwd(base_dir)


OUT_dir <- file.path(base_dir, "geomxtools")
QC_dir <- file.path(base_dir, "geomxtools", "QC")
plot_dir <- file.path(base_dir, "geomxtools", "results", "plots")

datadir <- system.file("extdata", "WTA_NGS_Example",
                       package="GeoMxWorkflows")

DCCFiles <- dir(file.path(base_dir, "raw_data/DCC-20250203"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles = file.path(base_dir, "ref/pkcs", "Hs_R_NGS_WTA_v1.0.pkc")

SampleAnnotationFile <- file.path(base_dir, "raw_data/annotations", "annotations.xlsx")

demoData <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = "Template",
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("Aoi", "Roi"),
                         experimentDataColNames = c("(v1.0) Human NGS Whole Transcriptome Atlas RNA"))

## 3.Study Design

library(knitr)
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

library(dplyr)
library(ggforce)
library(networkD3)

sankeyCols <- c("source", "target", "value")

link1 <- count(pData(demoData), individual, `slide name`)
link2 <- count(pData(demoData), `slide name`, class)
link3 <- count(pData(demoData),  class, region)
#link3 <- count(pData(demoData),  region, Segment)

colnames(link1) <- sankeyCols
colnames(link2) <- sankeyCols
colnames(link3) <- sankeyCols

links <- rbind(link1,link2,link3)
nodes <- unique(data.frame(name=c(links$source, links$target)))

# sankeyNetwork is 0 based, not 1 based
links$source <- as.integer(match(links$source,nodes$name)-1)
links$target <- as.integer(match(links$target,nodes$name)-1)


sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                        Target = "target", Value = "value", NodeID = "name",
                        units = "TWh", fontSize = 18, nodeWidth = 30)
# Save as an interactive HTML file
saveWidget(sankey, file.path(plot_dir, "sankey_plot.html"), selfcontained = TRUE)

webshot(file.path(plot_dir, "sankey_plot.html"), 
        file.path(plot_dir, "sankey_plot.png"))


## QC & Pre-processing

# Shift counts to one
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

bu = demoData

QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 100,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
demoData <-
  setSegmentQCFlags(demoData, 
                    qcCutoffs = QC_params)

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))


library(ggplot2)

col_by <- "Segment"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL,
                         save_path = NULL) {  # New argument for saving the plot
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if (!is.null(scale_trans)) {
    plt <- plt + scale_x_continuous(trans = scale_trans)
  }
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = plt, width = 8, height = 6, dpi = 300)
  }
  return(plt)  # Return the plot object for further use if needed
}


QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80,
             save_path = file.path(QC_dir, "QC.trimmed.png"))

QC_histogram(sData(demoData), "Aligned (%)", col_by, 80,
             save_path = file.path(QC_dir, "QC.aligned.png"))

QC_histogram(sData(demoData), "Saturated (%)", col_by, 50,
             save_path = file.path(QC_dir, "QC.saturated.png")) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")

QC_histogram(sData(demoData), "Area", col_by, 1000, scale_trans = "log10",
             save_path = file.path(QC_dir, "QC.area.png"))

QC_histogram(sData(demoData), "Nuclei", col_by, 100,
             save_path = file.path(QC_dir, "QC.nuclei.png"))


# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(demoData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}


# detatch neg_geomean columns ahead of aggregateCounts call
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(demoData)$NTC),
      col.names = c("NTC Count", "# of Segments"))
write.table(table(NTC_Count = sData(demoData)$NTC), 
            file.path(QC_dir, "NTC_Count.txt"),
            quote = F, sep = "\t", col.names = T, row.names = T)

kable(QC_Summary, caption = "QC Summary Table for each Segment")
write.table(table(NTC_Count = sData(demoData)$NTC), 
            file.path(QC_dir, "NTC_Count.txt"),
            quote = F, sep = "\t", col.names = T, row.names = T)


demoData <- demoData[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(demoData)
#Features  Samples 
#18815       91

### 4.2 Probe QC
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(demoData, 
         fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
#> Features  Samples 
#>    18815      91
demoData <- ProbeQCPassed


### 4.3 Create Gene-level Count Data
# Check how many unique targets the object has
length(unique(featureData(demoData)[["TargetName"]]))
#> [1] 18677

# collapse to targets
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)
#> Features  Samples 
#>    18677      91
exprs(target_demoData)[1:5, 1:2]


### 4.4 Limit of Quantification
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] * 
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}
pData(target_demoData)$LOQ <- LOQ


### 4.5 Filtering
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_demoData)$Module == module
  Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]


# Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
plt <- ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = region)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")
ggsave(filename = file.path(QC_dir, "Segment_Gene_Detection.png"), 
       plot = plt, width = 6, height = 6, dpi = 300)

# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_demoData)$DetectionThreshold,
            pData(target_demoData)$class))
write.table(table(pData(target_demoData)$DetectionThreshold,
                  pData(target_demoData)$class), 
            file.path(QC_dir, "Gene_Detection_Rate.txt"),
            quote = F, sep = "\t", col.names = T, row.names = T)

target_demoData <-
#  target_demoData[, pData(target_demoData)$GeneDetectionRate >= .1]
  target_demoData[, pData(target_demoData)$GeneDetectionRate >= .05]

dim(target_demoData)
#Features  Samples 
#18677       90

library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

# Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_demoData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))



# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
rownames(plot_detect) <- plot_detect$Freq

plt <- ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")
ggsave(filename = file.path(QC_dir, "Total_Detected_Genes.png"), 
       plot = plt, width = 6, height = 6, dpi = 300)

# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_demoData <- 
  target_demoData[fData(target_demoData)$DetectionRate >= 0.1 |
                    fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData)
#> Features  Samples 
#>    10131      221

# retain only detected genes of interest
goi <- goi[goi %in% rownames(target_demoData)]


### 5 Normalization
library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "region"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_demoData)),
             Segment = colnames(exprs(target_demoData)),
             Annotation = pData(target_demoData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_demoData), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_demoData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plt <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
ggsave(filename = file.path(QC_dir, "Normalization.png"), 
       plot = plt, width = 8, height = 6, dpi = 300)


# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_demoData <- normalize(target_demoData ,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_demoData <- normalize(target_demoData ,
                             norm_method = "neg", 
                             fromElt = "exprs",
                             toElt = "neg_norm")

assayDataElement(object = target_demoData, elt = "log_q") <-
  assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")


raw_counts <- assayDataElement(target_demoData, "exprs")  # assumes these are raw counts

dge <- DGEList(counts = raw_counts)
dge <- calcNormFactors(dge, method = "TMM")
pData(target_demoData)$TMM_factor <- dge$samples$norm.factors
tmm_cpm <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)  # linear scale
log_tmm_cpm <- cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
assayDataElement(target_demoData, "tmm_cpm") <- tmm_cpm
assayDataElement(target_demoData, "log_tmm_cpm") <- log_tmm_cpm


saveRDS(target_demoData, 
        file = file.path(OUT_dir, "target_demoData.norm.rds"))


# visualize the first 10 segments with each normalization method
png(file.path(QC_dir, "Raw_Counts.png"), width = 800, height = 600, res = 150)
boxplot(exprs(target_demoData)[,1:10],
         col = "#9EDAE5", main = "Raw Counts",
         log = "y", names = 1:10, xlab = "Segment",
         ylab = "Counts, Raw")
dev.off()

png(file.path(QC_dir, "Q3_Norm_Counts.png"), width = 800, height = 600, res = 150)
boxplot(assayDataElement(target_demoData[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")
dev.off()

png(file.path(QC_dir, "Neg_Norm_Counts.png"), width = 800, height = 600, res = 150)
boxplot(assayDataElement(target_demoData[,1:10], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Neg. Normalized")
dev.off()


### 6 Unsupervised Analysis
library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),  
       config = custom_umap)
#> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
#> Also defined by 'spam'
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
plt <- ggplot(pData(target_demoData),
       aes(x = UMAP1, y = UMAP2, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()
ggsave(filename = file.path(QC_dir, "UMAP.png"), 
       plot = plt, width = 6, height = 6, dpi = 300)

# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),
        perplexity = ncol(target_demoData)*.15)
pData(target_demoData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
plt <- ggplot(pData(target_demoData),
       aes(x = tSNE1, y = tSNE2, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()
ggsave(filename = file.path(QC_dir, "tSNE.png"), 
       plot = plt, width = 6, height = 6, dpi = 300)


library(pheatmap)  # for pheatmap
# create a log2 transform of the data for analysis
assayDataElement(object = target_demoData, elt = "log_q") <-
  assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_demoData,
                         elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]
#>   CAMK2N1    AKR1C1      AQP2     GDF15       REN 
#> 0.5886006 0.5114973 0.4607206 0.4196469 0.4193216

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
plt <- pheatmap(assayDataElement(target_demoData[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
           pData(target_demoData)[, c("class", "Segment", "region")])
ggsave(filename = file.path(QC_dir, "Clustering_CV_Genes.png"), 
       plot = plt, width = 6, height = 6, dpi = 300)

save.image(file=file.path(OUT_dir, 'geomxtools.RData'))
