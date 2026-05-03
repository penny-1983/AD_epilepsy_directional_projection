# =========================================================
# 00_setup.R
# Project setup
# =========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

proj_dir <- getwd()

dir_raw <- file.path(proj_dir, "data", "raw")
dir_signature <- file.path(proj_dir, "data", "directional_signature")
dir_results <- file.path(proj_dir, "results")
dir_tables <- file.path(dir_results, "source_tables")
dir_figdata <- file.path(dir_results, "figure_source_data")

dir.create(dir_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_signature, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_results, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figdata, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "edgeR",
  "limma",
  "GSVA",
  "singscore",
  "fgsea",
  "msigdbr",
  "ggplot2",
  "pheatmap",
  "STRINGdb",
  "igraph",
  "openxlsx"
)

message("Project directory: ", proj_dir)
message("Please make sure required packages are installed before running the scripts.")
