# =========================================================
# 02_project_signature_GSE5281.R
# Purpose:
# Prepare GSE5281 gene-level expression matrix and project
# the direction-aware epilepsy-derived signature onto AD tissues.
# =========================================================

source("scripts/00_setup.R")

# -------------------------
# Required packages
# -------------------------
suppressPackageStartupMessages({
  library(GSVA)
})

# -------------------------
# Expected input files
# -------------------------
# Users should prepare:
# 1. GSE5281 probe-level expression matrix:
#    data/raw/GSE5281/GSE5281_expression_probe_matrix.csv
#    - first column: Probe
#    - remaining columns: sample IDs
#
# 2. GSE5281 sample metadata:
#    data/raw/GSE5281/GSE5281_sample_metadata.csv
#    required columns: Sample, Diagnosis, BrainRegion
#
# 3. GPL570 annotation file:
#    data/raw/GSE5281/GPL570-55999.txt
#
# 4. Directional signatures:
#    data/directional_signature/Epilepsy_Up_top100.txt
#    data/directional_signature/Epilepsy_Down_top100.txt

expr_file <- file.path(dir_raw, "GSE5281", "GSE5281_expression_probe_matrix.csv")
meta_file <- file.path(dir_raw, "GSE5281", "GSE5281_sample_metadata.csv")
annot_file <- file.path(dir_raw, "GSE5281", "GPL570-55999.txt")

up_file <- file.path(dir_signature, "Epilepsy_Up_top100.txt")
down_file <- file.path(dir_signature, "Epilepsy_Down_top100.txt")

required_files <- c(expr_file, meta_file, annot_file, up_file, down_file)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

# -------------------------
# Step 1. Read expression and metadata
# -------------------------
expr_df <- read.csv(
  expr_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- read.csv(
  meta_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_meta_cols <- c("Sample", "Diagnosis", "BrainRegion")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    "Metadata file must contain columns: ",
    paste(required_meta_cols, collapse = ", ")
  )
}

meta$Sample <- trimws(as.character(meta$Sample))
meta$Diagnosis <- trimws(as.character(meta$Diagnosis))
meta$BrainRegion <- trimws(as.character(meta$BrainRegion))

# -------------------------
# Step 2. Prepare probe-level matrix
# -------------------------
if (!"Probe" %in% colnames(expr_df)) {
  stop("Expression matrix must contain a first column named 'Probe'.")
}

rownames(expr_df) <- toupper(trimws(expr_df$Probe))
expr_mat <- as.matrix(expr_df[, setdiff(colnames(expr_df), "Probe"), drop = FALSE])
mode(expr_mat) <- "numeric"

# -------------------------
# Step 3. Probe-to-gene mapping using GPL570 annotation
# -------------------------
annot <- read.delim(
  annot_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  skip = 16
)

if (!all(c("ID", "Gene Symbol") %in% colnames(annot))) {
  stop("GPL570 annotation file must contain columns 'ID' and 'Gene Symbol'.")
}

probe2symbol <- annot[, c("ID", "Gene Symbol")]
colnames(probe2symbol) <- c("Probe", "GeneSymbol")

probe2symbol <- probe2symbol[
  !is.na(probe2symbol$GeneSymbol) & probe2symbol$GeneSymbol != "",
]

probe2symbol$Probe <- toupper(trimws(probe2symbol$Probe))
probe2symbol$GeneSymbol <- sapply(
  strsplit(probe2symbol$GeneSymbol, " /// "),
  `[`,
  1
)
probe2symbol$GeneSymbol <- toupper(trimws(probe2symbol$GeneSymbol))

probe2symbol <- probe2symbol[
  !is.na(probe2symbol$GeneSymbol) & probe2symbol$GeneSymbol != "",
]

expr_probe <- data.frame(
  Probe = rownames(expr_mat),
  expr_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_annot <- merge(probe2symbol, expr_probe, by = "Probe")

sample_cols <- colnames(expr_mat)

expr_gene_df <- aggregate(
  x = expr_annot[, sample_cols, drop = FALSE],
  by = list(GeneSymbol = expr_annot$GeneSymbol),
  FUN = mean
)

rownames(expr_gene_df) <- expr_gene_df$GeneSymbol
expr_gene <- as.matrix(expr_gene_df[, -1, drop = FALSE])
mode(expr_gene) <- "numeric"

write.csv(
  data.frame(Gene = rownames(expr_gene), expr_gene, check.names = FALSE),
  file = file.path(dir_tables, "GSE5281_gene_level_expression_matrix.csv"),
  row.names = FALSE
)

message("Gene-level expression matrix generated: ",
        nrow(expr_gene), " genes x ", ncol(expr_gene), " samples.")

# -------------------------
# Step 4. Read directional signature
# -------------------------
read_signature <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  x <- x[x != "" & !is.na(x)]
  unique(toupper(x))
}

up_sig <- read_signature(up_file)
down_sig <- read_signature(down_file)

matched_up <- intersect(up_sig, rownames(expr_gene))
matched_down <- intersect(down_sig, rownames(expr_gene))

message("Matched SOZ-up genes: ", length(matched_up))
message("Matched SOZ-down genes: ", length(matched_down))

if (length(matched_up) < 10 || length(matched_down) < 10) {
  stop("Too few matched signature genes. Please check gene symbols.")
}

# -------------------------
# Step 5. ssGSEA projection
# -------------------------
gene_sets <- list(
  Up = matched_up,
  Down = matched_down
)

ssgsea_res <- gsva(
  expr = expr_gene,
  gset.idx.list = gene_sets,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE,
  verbose = FALSE
)

score_df <- data.frame(
  Sample = colnames(expr_gene),
  ssGSEA_Up = as.numeric(ssgsea_res["Up", ]),
  ssGSEA_Down = as.numeric(ssgsea_res["Down", ]),
  stringsAsFactors = FALSE
)

score_df$EpilepsyScore <- score_df$ssGSEA_Up - score_df$ssGSEA_Down
score_df$Sample <- trimws(as.character(score_df$Sample))

score_df2 <- merge(score_df, meta, by = "Sample")

score_df2$Diagnosis <- factor(score_df2$Diagnosis, levels = c("Control", "AD"))
score_df2$BrainRegion <- factor(score_df2$BrainRegion)

write.csv(
  score_df2,
  file = file.path(dir_tables, "GSE5281_EpilepsyScore_top100.csv"),
  row.names = FALSE
)

# -------------------------
# Step 6. Overall AD vs Control comparison
# -------------------------
overall_wilcox <- wilcox.test(EpilepsyScore ~ Diagnosis, data = score_df2)
overall_ttest <- t.test(EpilepsyScore ~ Diagnosis, data = score_df2)

overall_summary <- aggregate(
  EpilepsyScore ~ Diagnosis,
  data = score_df2,
  function(x) {
    c(
      n = length(x),
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE)
    )
  }
)

overall_summary <- do.call(data.frame, overall_summary)

overall_tests <- data.frame(
  Wilcoxon_P = overall_wilcox$p.value,
  Welch_t_P = overall_ttest$p.value
)

write.csv(
  overall_summary,
  file = file.path(dir_tables, "GSE5281_overall_EpilepsyScore_summary.csv"),
  row.names = FALSE
)

write.csv(
  overall_tests,
  file = file.path(dir_tables, "GSE5281_overall_EpilepsyScore_tests.csv"),
  row.names = FALSE
)

# -------------------------
# Step 7. Region-wise comparisons
# -------------------------
region_list <- by(score_df2, score_df2$BrainRegion, function(df) {
  if (length(unique(df$Diagnosis)) < 2) {
    return(NULL)
  }
  
  wt <- wilcox.test(EpilepsyScore ~ Diagnosis, data = df)
  tt <- t.test(EpilepsyScore ~ Diagnosis, data = df)
  
  ad_vals <- df$EpilepsyScore[df$Diagnosis == "AD"]
  ctrl_vals <- df$EpilepsyScore[df$Diagnosis == "Control"]
  
  data.frame(
    Region = unique(df$BrainRegion),
    AD_n = length(ad_vals),
    Control_n = length(ctrl_vals),
    AD_mean = mean(ad_vals, na.rm = TRUE),
    Control_mean = mean(ctrl_vals, na.rm = TRUE),
    AD_median = median(ad_vals, na.rm = TRUE),
    Control_median = median(ctrl_vals, na.rm = TRUE),
    Wilcoxon_P = wt$p.value,
    Welch_t_P = tt$p.value,
    Direction = ifelse(
      mean(ad_vals, na.rm = TRUE) < mean(ctrl_vals, na.rm = TRUE),
      "AD < Control",
      "AD > Control"
    ),
    stringsAsFactors = FALSE
  )
})

region_stats <- do.call(rbind, region_list)
region_stats$FDR_Wilcoxon <- p.adjust(region_stats$Wilcoxon_P, method = "BH")
region_stats <- region_stats[order(region_stats$Wilcoxon_P), ]

write.csv(
  region_stats,
  file = file.path(dir_tables, "GSE5281_region_wise_EpilepsyScore_comparison.csv"),
  row.names = FALSE
)

# -------------------------
# Step 8. Region-adjusted linear model
# -------------------------
fit_main <- lm(EpilepsyScore ~ Diagnosis + BrainRegion, data = score_df2)
fit_int <- lm(EpilepsyScore ~ Diagnosis * BrainRegion, data = score_df2)

model_summary <- capture.output({
  cat("Main model: EpilepsyScore ~ Diagnosis + BrainRegion\n")
  print(summary(fit_main))
  cat("\nInteraction model: EpilepsyScore ~ Diagnosis * BrainRegion\n")
  print(summary(fit_int))
  cat("\nModel comparison:\n")
  print(anova(fit_main, fit_int))
})

writeLines(
  model_summary,
  con = file.path(dir_tables, "GSE5281_linear_model_summary.txt")
)

message("GSE5281 projection and region-wise analysis completed.")
