library(DESeq2)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)

if (length(script_arg) > 0) {
  script_path <- sub("^--file=", "", script_arg[1])
  setwd(normalizePath(file.path(dirname(script_path), "..", "..")))
}

input_directory <- Sys.getenv("FEATURECOUNTS_DIR", unset = file.path("datos", "featurecounts"))
if (!dir.exists(input_directory)) {
  stop("Configure FEATURECOUNTS_DIR con la carpeta de matrices FeatureCounts.")
}

output_directory <- file.path("resultados", "DEGs", "deseq2")
graphics_directory <- file.path("resultados", "graficos", "deseq2")

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(graphics_directory, recursive = TRUE, showWarnings = FALSE)

# Variables de análisis
select_regions <- c("Magallanes") # Quitar según que región se use(Aysen, Magallanes)
select_seasons <- c("WIN", "SUM")  # Agregar o quitar estaciones del año (AUT, WIN, SPR, SUM)

# Cambia según tipo de análisis (season si hay más de una estación o region si hay más de una región)
if (length(select_seasons) == 2 & length(select_regions) == 1) {
  analysis_type <- "season"
} else if (length(select_regions) == 2 & length(select_seasons) == 1) {
  analysis_type <- "region"
} else {
  stop("Aplique bien el formato: 2 regiones o 2 estaciones.")
}

aysen <- read.table(file.path(input_directory, "counts_Aysen.txt"), header=TRUE, row.names=1, check.names=FALSE)
magallanes <- read.table(file.path(input_directory, "counts_Magallanes.txt"), header=TRUE, row.names=1, check.names=FALSE)

aysen_clean <- aysen[, 6:ncol(aysen)]
magallanes_clean <- magallanes[, 6:ncol(magallanes)]

colnames(aysen_clean) <- c(
  "P102_A_SPR",
  "P10_A_SUM",
  "P1110_A_SPR",
  "P11_A_SUM",
  "P1213_A_SPR",
  "P29_A_SUM",
  "P41_A_WIN",
  "P56_A_WIN",
  "P611_A_WIN"
)

colnames(magallanes_clean) <- c(
  "P11_M_SUM",
  "P1_M_AUT",
  "P311_M_SUM",
  "P3_M_AUT",
  "P4_M_WIN",
  "P5_M_WIN",
  "P6_M_WIN",
  "P7_M_SPR",
  "P8_M_SPR",
  "P9_M_SPR"
)

counts_matrix <- merge(aysen_clean, magallanes_clean, by="row.names")
rownames(counts_matrix) <- counts_matrix$Row.names

counts_matrix$Row.names <- NULL

samples <- colnames(counts_matrix)

region <- ifelse(
  grepl("_A_", samples),
  "Aysen",
  "Magallanes"
)

season <- sub(".*_", "", samples)

metadata <- data.frame(
  sample = samples,
  region = region,
  season = season
)

rownames(metadata) <- metadata$sample

metadata

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData = metadata,
  design = ~ region + season
)

dds_filter <- dds[
  ,
  dds$region %in% select_regions & #Filtrar según región
  dds$season %in% select_seasons   #Filtrar según estación 
]

dds_filter$season <- droplevels(dds_filter$season)
dds_filter$region <- droplevels(dds_filter$region)

#Según el tipo de analisis aplicamos distintos contrastes, y cambiamos el título para el gráfico y csv.
if (analysis_type == "season") {
  design(dds_filter) <- ~ season
  aux_contrast <- c(
    "season",
    select_seasons[1],
    select_seasons[2]
  )

  title_graphic <- paste(
    select_regions[1],
    ":",
    select_seasons[1],
    "vs",
    select_seasons[2]
  )

  name_csv <- paste0(
    "DEGs_",
    select_regions[1],
    "_",
    select_seasons[1],
    "_vs_",
    select_seasons[2],
    ".csv"
  )

} else if (analysis_type == "region") {
  design(dds_filter) <- ~ region

  aux_contrast <- c(
    "region",
    select_regions[1],
    select_regions[2]
  )

  title_graphic <- paste(
    select_regions[1],
    "vs",
    select_regions[2],
    "-",
    select_seasons[1]
  )

  name_csv <- paste0(
    "DEGs_",
    select_regions[1],
    "_vs_",
    select_regions[2],
    "_",
    select_seasons[1],
    ".csv"
  )
}

keep <- rowSums(counts(dds_filter) >= 10) >= 2

dds_filter <- dds_filter[keep, ]
dds_filter <- DESeq(dds_filter)

res <- results(dds_filter, alpha=0.05, contrast = aux_contrast)
res <- res[!is.na(res$padj), ]

res_filter <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 2, ]

write.csv(as.data.frame(res_filter), file.path(output_directory, name_csv))

library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = title_graphic,
                pointSize = 2.0,
                labSize = 3.0,
                col = c('black', 'black', 'black', 'red'),
                xlim =c(-13,13),
                ylim = c(0,161)
                )

vsd <- vst(dds_filter, blind = FALSE)

plotPCA(
  vsd,
  intgroup = ifelse(
    analysis_type == "season",
    "season",
    "region"
  )
)

ggsave(file.path(graphics_directory, "PCA_A_SPR_SUM.png"), width = 6, height = 3, units = "in", dpi = 300)

library(pheatmap)

sig_res <- res[
  res$padj < 0.05 &
    abs(res$log2FoldChange) > 2,
]

top_sig <- head(
  sig_res[order(sig_res$padj), ],
  30
)

sig_genes <- rownames(top_sig)

heatmap_matrix <- assay(vsd)[sig_genes, ]

pheatmap(
  heatmap_matrix,
  scale = "row",
  show_rownames = TRUE,
  annotation_col = as.data.frame(colData(vsd)[, c("region", "season")]),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = paste("Heatmap -", title_graphic)
)
