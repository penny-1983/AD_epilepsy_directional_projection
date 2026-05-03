# =========================================================
# 05_external_cohort_evaluation.R
# Purpose:
# External evaluation of the direction-aware epilepsy-derived
# projection score in independent AD cohorts:
#   - GSE132903: MTG-focused cohort
#   - GSE48350: multi-region cohort
# =========================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(GSVA)
  library(ggplot2)
})

# lme4 is optional and only required for mixed-effects sensitivity analysis
has_lme4 <- requireNamespace("lme4", quietly = TRUE)

# -------------------------
# Helper functions
# -------------------------
read_signature <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  x <- x[x != "" & !is.na(x)]
  unique(toupper(x))
}

read_gene_matrix <- function(file) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  
  if (!"Gene" %in% colnames(df)) {
    stop("Expression matrix must contain a first column named 'Gene': ", file)
  }
  
  mat <- as.matrix(df[, setdiff(colnames(df), "Gene"), drop = FALSE])
  rownames(mat) <- toupper(trimws(df$Gene))
  mode(mat) <- "numeric"
  
  mat
}

read_metadata <- function(file) {
  meta <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  
  required_cols <- c("Sample", "Diagnosis")
  miss <- setdiff(required_cols, colnames(meta))
  if (length(miss) > 0) {
    stop("Metadata file is missing required columns: ", paste(miss, collapse = ", "))
  }
  
  meta$Sample <- trimws(as.character(meta$Sample))
  meta$Diagnosis <- trimws(as.character(meta$Diagnosis))
  
  if ("BrainRegion" %in% colnames(meta)) {
    meta$BrainRegion <- trimws(as.character(meta$BrainRegion))
  } else {
    meta$BrainRegion <- "Unknown"
  }
  
  if ("SubjectID" %in% colnames(meta)) {
    meta$SubjectID <- trimws(as.character(meta$SubjectID))
  }
  
  if ("Age" %in% colnames(meta)) {
    meta$Age <- suppressWarnings(as.numeric(meta$Age))
  }
  
  if ("Sex" %in% colnames(meta)) {
    meta$Sex <- factor(trimws(as.character(meta$Sex)))
  }
  
  meta$Diagnosis <- factor(meta$Diagnosis, levels = c("Control", "AD"))
  meta$BrainRegion <- factor(meta$BrainRegion)
  
  meta
}

calc_projection_score <- function(expr_mat, up_genes, down_genes) {
  up_use <- intersect(up_genes, rownames(expr_mat))
  down_use <- intersect(down_genes, rownames(expr_mat))
  
  if (length(up_use) < 10 || length(down_use) < 10) {
    stop("Too few matched genes for projection.")
  }
  
  gene_sets <- list(
    Up = up_use,
    Down = down_use
  )
  
  ss <- gsva(
    expr = expr_mat,
    gset.idx.list = gene_sets,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = TRUE,
    verbose = FALSE
  )
  
  score <- as.numeric(ss["Up", ] - ss["Down", ])
  names(score) <- colnames(expr_mat)
  
  list(
    score = score,
    matched_up = up_use,
    matched_down = down_use
  )
}

merge_score_metadata <- function(expr_mat, meta, score) {
  df <- data.frame(
    Sample = names(score),
    EpilepsyScore = as.numeric(score),
    stringsAsFactors = FALSE
  )
  
  df <- merge(df, meta, by = "Sample")
  df <- df[df$Diagnosis %in% c("Control", "AD"), , drop = FALSE]
  df$Diagnosis <- factor(as.character(df$Diagnosis), levels = c("Control", "AD"))
  
  df
}

overall_stats <- function(df, cohort_name, subset_name = "Overall") {
  if (length(unique(df$Diagnosis)) < 2) return(NULL)
  
  wt <- wilcox.test(EpilepsyScore ~ Diagnosis, data = df)
  tt <- t.test(EpilepsyScore ~ Diagnosis, data = df)
  
  ad_vals <- df$EpilepsyScore[df$Diagnosis == "AD"]
  ctrl_vals <- df$EpilepsyScore[df$Diagnosis == "Control"]
  
  data.frame(
    Cohort = cohort_name,
    Subset = subset_name,
    N_total = nrow(df),
    N_Control = length(ctrl_vals),
    N_AD = length(ad_vals),
    Mean_Control = mean(ctrl_vals, na.rm = TRUE),
    Mean_AD = mean(ad_vals, na.rm = TRUE),
    Median_Control = median(ctrl_vals, na.rm = TRUE),
    Median_AD = median(ad_vals, na.rm = TRUE),
    Direction = ifelse(mean(ad_vals, na.rm = TRUE) < mean(ctrl_vals, na.rm = TRUE),
                       "AD < Control", "AD > Control"),
    Wilcoxon_P = wt$p.value,
    Welch_t_P = tt$p.value,
    stringsAsFactors = FALSE
  )
}

region_stats <- function(df, cohort_name) {
  regions <- unique(as.character(df$BrainRegion))
  
  out <- lapply(regions, function(rg) {
    sub <- df[df$BrainRegion == rg, , drop = FALSE]
    if (length(unique(sub$Diagnosis)) < 2) return(NULL)
    
    res <- overall_stats(sub, cohort_name, subset_name = rg)
    res
  })
  
  out <- do.call(rbind, out)
  if (!is.null(out)) {
    out$FDR_Wilcoxon <- p.adjust(out$Wilcoxon_P, method = "BH")
  }
  
  out
}

run_region_adjusted_lm <- function(df, cohort_name, model_name = "Region-adjusted linear model") {
  if (!"BrainRegion" %in% colnames(df)) return(NULL)
  if (length(unique(df$Diagnosis)) < 2) return(NULL)
  
  fit <- lm(EpilepsyScore ~ Diagnosis + BrainRegion, data = df)
  
  if (!"DiagnosisAD" %in% rownames(coef(summary(fit)))) return(NULL)
  
  beta <- coef(summary(fit))["DiagnosisAD", "Estimate"]
  se <- coef(summary(fit))["DiagnosisAD", "Std. Error"]
  p <- coef(summary(fit))["DiagnosisAD", "Pr(>|t|)"]
  
  data.frame(
    Cohort = cohort_name,
    Model = model_name,
    Diagnosis_beta = beta,
    CI_lower = beta - 1.96 * se,
    CI_upper = beta + 1.96 * se,
    P_value = p,
    Direction = ifelse(beta < 0, "AD < Control", "AD > Control"),
    stringsAsFactors = FALSE
  )
}

run_mixed_model <- function(df, cohort_name, adjusted = FALSE) {
  if (!has_lme4) {
    message("Package lme4 is not installed; mixed-effects analysis will be skipped.")
    return(NULL)
  }
  
  if (!"SubjectID" %in% colnames(df)) {
    message("SubjectID not found; mixed-effects analysis will be skipped.")
    return(NULL)
  }
  
  df <- df[!is.na(df$SubjectID), , drop = FALSE]
  if (length(unique(df$SubjectID)) < 5) return(NULL)
  
  if (adjusted && all(c("Age", "Sex") %in% colnames(df))) {
    df2 <- df[complete.cases(df[, c("EpilepsyScore", "Diagnosis", "BrainRegion", "Age", "Sex", "SubjectID")]), ]
    model_formula <- EpilepsyScore ~ Diagnosis + BrainRegion + Age + Sex + (1 | SubjectID)
    model_name <- "Mixed-effects model adjusted for brain region, age, and sex"
  } else {
    df2 <- df[complete.cases(df[, c("EpilepsyScore", "Diagnosis", "BrainRegion", "SubjectID")]), ]
    model_formula <- EpilepsyScore ~ Diagnosis + BrainRegion + (1 | SubjectID)
    model_name <- "Mixed-effects model with subject-specific random intercept"
  }
  
  if (nrow(df2) < 10 || length(unique(df2$Diagnosis)) < 2) return(NULL)
  
  fit <- lme4::lmer(model_formula, data = df2)
  sm <- summary(fit)$coefficients
  
  if (!"DiagnosisAD" %in% rownames(sm)) return(NULL)
  
  beta <- sm["DiagnosisAD", "Estimate"]
  se <- sm["DiagnosisAD", "Std. Error"]
  tval <- sm["DiagnosisAD", "t value"]
  
  data.frame(
    Cohort = cohort_name,
    Model = model_name,
    N = nrow(df2),
    N_subjects = length(unique(df2$SubjectID)),
    Diagnosis_beta = beta,
    CI_lower = beta - 1.96 * se,
    CI_upper = beta + 1.96 * se,
    Approx_t = tval,
    Direction = ifelse(beta < 0, "AD < Control", "AD > Control"),
    stringsAsFactors = FALSE
  )
}

run_subject_collapsed <- function(df, cohort_name) {
  if (!"SubjectID" %in% colnames(df)) {
    message("SubjectID not found; subject-collapsed analysis will be skipped.")
    return(NULL)
  }
  
  df <- df[!is.na(df$SubjectID), , drop = FALSE]
  
  collapsed <- aggregate(
    EpilepsyScore ~ SubjectID + Diagnosis,
    data = df,
    FUN = mean
  )
  
  if (length(unique(collapsed$Diagnosis)) < 2) return(NULL)
  
  wt <- wilcox.test(EpilepsyScore ~ Diagnosis, data = collapsed)
  tt <- t.test(EpilepsyScore ~ Diagnosis, data = collapsed)
  fit <- lm(EpilepsyScore ~ Diagnosis, data = collapsed)
  
  beta <- coef(summary(fit))["DiagnosisAD", "Estimate"]
  se <- coef(summary(fit))["DiagnosisAD", "Std. Error"]
  p <- coef(summary(fit))["DiagnosisAD", "Pr(>|t|)"]
  
  data.frame(
    Cohort = cohort_name,
    Model = "Subject-collapsed sensitivity analysis",
    N_subjects = nrow(collapsed),
    Diagnosis_beta = beta,
    CI_lower = beta - 1.96 * se,
    CI_upper = beta + 1.96 * se,
    P_value = p,
    Wilcoxon_P = wt$p.value,
    Welch_t_P = tt$p.value,
    Direction = ifelse(beta < 0, "AD < Control", "AD > Control"),
    stringsAsFactors = FALSE
  )
}

plot_external_boxplot <- function(df, out_file, title_text) {
  p <- ggplot(df, aes(x = Diagnosis, y = EpilepsyScore, fill = Diagnosis)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
    scale_fill_manual(values = c("Control" = "#F8766D", "AD" = "#00BFC4")) +
    labs(
      title = title_text,
      x = "Diagnosis",
      y = "EpilepsyScore"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  ggsave(out_file, p, width = 5.5, height = 4.5, dpi = 300)
}

# =========================================================
# Input signature
# =========================================================

up_file <- file.path(dir_signature, "Epilepsy_Up_top100.txt")
down_file <- file.path(dir_signature, "Epilepsy_Down_top100.txt")

if (!file.exists(up_file) || !file.exists(down_file)) {
  stop("Top100 directional signature files are missing.")
}

up_sig <- read_signature(up_file)
down_sig <- read_signature(down_file)

# =========================================================
# Part A. GSE132903 MTG-focused external evaluation
# =========================================================

gse132903_expr_file <- file.path(dir_raw, "GSE132903", "GSE132903_gene_level_expression_matrix.csv")
gse132903_meta_file <- file.path(dir_raw, "GSE132903", "GSE132903_sample_metadata.csv")

if (file.exists(gse132903_expr_file) && file.exists(gse132903_meta_file)) {
  
  expr_132903 <- read_gene_matrix(gse132903_expr_file)
  meta_132903 <- read_metadata(gse132903_meta_file)
  
  proj_132903 <- calc_projection_score(expr_132903, up_sig, down_sig)
  df_132903 <- merge_score_metadata(expr_132903, meta_132903, proj_132903$score)
  
  write.csv(
    df_132903,
    file = file.path(dir_tables, "GSE132903_EpilepsyScore_external.csv"),
    row.names = FALSE
  )
  
  stats_132903 <- overall_stats(df_132903, "GSE132903", "MTG-focused cohort")
  
  write.csv(
    stats_132903,
    file = file.path(dir_tables, "GSE132903_external_MTG_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Cohort = "GSE132903",
      Matched_SOZ_up_genes = length(proj_132903$matched_up),
      Matched_SOZ_down_genes = length(proj_132903$matched_down)
    ),
    file = file.path(dir_tables, "GSE132903_matched_signature_genes.csv"),
    row.names = FALSE
  )
  
  write.csv(
    df_132903,
    file = file.path(dir_figdata, "Supplementary_Figure_S4_source_data.csv"),
    row.names = FALSE
  )
  
  plot_external_boxplot(
    df_132903,
    out_file = file.path(dir_results, "Supplementary_Figure_S4_GSE132903_external_MTG.png"),
    title_text = "GSE132903 external MTG evaluation"
  )
  
  message("GSE132903 external evaluation completed.")
  
} else {
  message("GSE132903 input files not found. Skipping GSE132903 external evaluation.")
}

# =========================================================
# Part B. GSE48350 multi-region external evaluation
# =========================================================

gse48350_expr_file <- file.path(dir_raw, "GSE48350", "GSE48350_gene_level_expression_matrix.csv")
gse48350_meta_file <- file.path(dir_raw, "GSE48350", "GSE48350_sample_metadata.csv")

if (file.exists(gse48350_expr_file) && file.exists(gse48350_meta_file)) {
  
  expr_48350 <- read_gene_matrix(gse48350_expr_file)
  meta_48350 <- read_metadata(gse48350_meta_file)
  
  proj_48350 <- calc_projection_score(expr_48350, up_sig, down_sig)
  df_48350 <- merge_score_metadata(expr_48350, meta_48350, proj_48350$score)
  
  write.csv(
    df_48350,
    file = file.path(dir_tables, "GSE48350_EpilepsyScore_external.csv"),
    row.names = FALSE
  )
  
  # Overall and region-wise summaries
  stats_48350_overall <- overall_stats(df_48350, "GSE48350", "Overall")
  stats_48350_region <- region_stats(df_48350, "GSE48350")
  lm_48350 <- run_region_adjusted_lm(df_48350, "GSE48350")
  
  write.csv(
    stats_48350_overall,
    file = file.path(dir_tables, "GSE48350_external_overall_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    stats_48350_region,
    file = file.path(dir_tables, "GSE48350_external_region_wise_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    lm_48350,
    file = file.path(dir_tables, "GSE48350_region_adjusted_linear_model.csv"),
    row.names = FALSE
  )
  
  # Age-restricted analysis if Age is available
  if ("Age" %in% colnames(df_48350)) {
    df_48350_age60 <- df_48350[!is.na(df_48350$Age) & df_48350$Age >= 60, ]
    
    if (nrow(df_48350_age60) > 0 && length(unique(df_48350_age60$Diagnosis)) == 2) {
      age60_overall <- overall_stats(df_48350_age60, "GSE48350", "Age >= 60")
      age60_lm <- run_region_adjusted_lm(
        df_48350_age60,
        "GSE48350",
        model_name = "Age >= 60 region-adjusted linear model"
      )
      
      write.csv(
        age60_overall,
        file = file.path(dir_tables, "GSE48350_age60_overall_summary.csv"),
        row.names = FALSE
      )
      
      write.csv(
        age60_lm,
        file = file.path(dir_tables, "GSE48350_age60_region_adjusted_linear_model.csv"),
        row.names = FALSE
      )
    }
  }
  
  # Dependence-aware mixed-effects sensitivity
  mixed_1 <- run_mixed_model(df_48350, "GSE48350", adjusted = FALSE)
  mixed_2 <- run_mixed_model(df_48350, "GSE48350", adjusted = TRUE)
  collapsed <- run_subject_collapsed(df_48350, "GSE48350")
  
  sensitivity_table <- do.call(
    rbind,
    Filter(Negate(is.null), list(lm_48350, mixed_1, mixed_2, collapsed))
  )
  
  if (!is.null(sensitivity_table)) {
    write.csv(
      sensitivity_table,
      file = file.path(dir_tables, "GSE48350_dependence_aware_sensitivity.csv"),
      row.names = FALSE
    )
    
    write.csv(
      sensitivity_table,
      file = file.path(dir_figdata, "Supplementary_Figure_S5_source_data.csv"),
      row.names = FALSE
    )
  }
  
  # Plot overall and region-wise
  plot_external_boxplot(
    df_48350,
    out_file = file.path(dir_results, "Figure6A_GSE48350_external_overall.png"),
    title_text = "GSE48350 external evaluation"
  )
  
  p_region <- ggplot(df_48350, aes(x = Diagnosis, y = EpilepsyScore, fill = Diagnosis)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.12, size = 1.4, alpha = 0.65) +
    facet_wrap(~ BrainRegion, scales = "free_y") +
    scale_fill_manual(values = c("Control" = "#F8766D", "AD" = "#00BFC4")) +
    labs(
      x = "Diagnosis",
      y = "EpilepsyScore"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(dir_results, "Figure6B_GSE48350_external_region_wise.png"),
    plot = p_region,
    width = 8.5,
    height = 5.5,
    dpi = 300
  )
  
  message("GSE48350 external evaluation completed.")
  
} else {
  message("GSE48350 input files not found. Skipping GSE48350 external evaluation.")
}

message("External cohort evaluation script completed.")
