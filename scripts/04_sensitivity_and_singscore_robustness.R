# =========================================================
# 04_sensitivity_and_singscore_robustness.R
# Purpose:
# 1) Evaluate signature-size robustness using top50/top100/top200 signatures.
# 2) Assess concordance between ssGSEA-based EpilepsyScore and centered singscore.
# =========================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(GSVA)
  library(singscore)
  library(ggplot2)
})

# -------------------------
# Input files
# -------------------------
expr_file <- file.path(dir_tables, "GSE5281_gene_level_expression_matrix.csv")
score_file <- file.path(dir_tables, "GSE5281_EpilepsyScore_top100.csv")

if (!file.exists(expr_file)) {
  stop("Missing gene-level expression matrix. Please run 02_project_signature_GSE5281.R first.")
}

if (!file.exists(score_file)) {
  stop("Missing top100 EpilepsyScore table. Please run 02_project_signature_GSE5281.R first.")
}

expr_df <- read.csv(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
score_df <- read.csv(score_file, check.names = FALSE, stringsAsFactors = FALSE)

if (!"Gene" %in% colnames(expr_df)) {
  stop("Expression matrix must contain a first column named 'Gene'.")
}

expr_gene <- as.matrix(expr_df[, setdiff(colnames(expr_df), "Gene"), drop = FALSE])
rownames(expr_gene) <- toupper(trimws(expr_df$Gene))
mode(expr_gene) <- "numeric"

score_df$Sample <- trimws(as.character(score_df$Sample))
score_df$Diagnosis <- factor(score_df$Diagnosis, levels = c("Control", "AD"))
score_df$BrainRegion <- factor(score_df$BrainRegion)

# -------------------------
# Helper functions
# -------------------------
read_signature <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  x <- x[x != "" & !is.na(x)]
  unique(toupper(x))
}

calc_ssgsea_score <- function(expr_mat, up_genes, down_genes) {
  up_use <- intersect(up_genes, rownames(expr_mat))
  down_use <- intersect(down_genes, rownames(expr_mat))
  
  if (length(up_use) < 10 || length(down_use) < 10) {
    stop("Too few matched genes for projection.")
  }
  
  gs <- list(Up = up_use, Down = down_use)
  
  ss <- gsva(
    expr = expr_mat,
    gset.idx.list = gs,
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

# =========================================================
# Part A. Signature-size sensitivity analysis
# =========================================================

signature_sizes <- c(50, 100, 200)

robust_list <- list()

for (n in signature_sizes) {
  up_file <- file.path(dir_signature, paste0("Epilepsy_Up_top", n, ".txt"))
  down_file <- file.path(dir_signature, paste0("Epilepsy_Down_top", n, ".txt"))
  
  if (!file.exists(up_file) || !file.exists(down_file)) {
    stop("Missing signature files for top", n)
  }
  
  up_sig <- read_signature(up_file)
  down_sig <- read_signature(down_file)
  
  proj <- calc_ssgsea_score(expr_gene, up_sig, down_sig)
  
  df <- score_df[, c("Sample", "Diagnosis", "BrainRegion")]
  df$EpilepsyScore <- as.numeric(proj$score[df$Sample])
  df$Signature <- paste0("Top", n)
  
  robust_list[[paste0("Top", n)]] <- list(
    data = df,
    matched_up = length(proj$matched_up),
    matched_down = length(proj$matched_down)
  )
}

# -------------------------
# Overall sensitivity table
# -------------------------
overall_sensitivity <- do.call(rbind, lapply(names(robust_list), function(nm) {
  df <- robust_list[[nm]]$data
  
  wt <- wilcox.test(EpilepsyScore ~ Diagnosis, data = df)
  tt <- t.test(EpilepsyScore ~ Diagnosis, data = df)
  fit <- lm(EpilepsyScore ~ Diagnosis + BrainRegion, data = df)
  
  data.frame(
    Signature = nm,
    Matched_SOZ_up_genes = robust_list[[nm]]$matched_up,
    Matched_SOZ_down_genes = robust_list[[nm]]$matched_down,
    Mean_Control = mean(df$EpilepsyScore[df$Diagnosis == "Control"], na.rm = TRUE),
    Mean_AD = mean(df$EpilepsyScore[df$Diagnosis == "AD"], na.rm = TRUE),
    Median_Control = median(df$EpilepsyScore[df$Diagnosis == "Control"], na.rm = TRUE),
    Median_AD = median(df$EpilepsyScore[df$Diagnosis == "AD"], na.rm = TRUE),
    Wilcoxon_P = wt$p.value,
    Welch_t_P = tt$p.value,
    Diagnosis_beta_adjusted_for_region = coef(summary(fit))["DiagnosisAD", "Estimate"],
    Diagnosis_beta_P = coef(summary(fit))["DiagnosisAD", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}))

write.csv(
  overall_sensitivity,
  file = file.path(dir_tables, "GSE5281_signature_size_overall_sensitivity.csv"),
  row.names = FALSE
)

# -------------------------
# Region-wise sensitivity table
# -------------------------
region_sensitivity <- do.call(rbind, lapply(names(robust_list), function(nm) {
  df <- robust_list[[nm]]$data
  
  out <- by(df, df$BrainRegion, function(sub) {
    if (length(unique(sub$Diagnosis)) < 2) {
      return(NULL)
    }
    
    wt <- wilcox.test(EpilepsyScore ~ Diagnosis, data = sub)
    tt <- t.test(EpilepsyScore ~ Diagnosis, data = sub)
    
    data.frame(
      Signature = nm,
      BrainRegion = unique(as.character(sub$BrainRegion)),
      Mean_Control = mean(sub$EpilepsyScore[sub$Diagnosis == "Control"], na.rm = TRUE),
      Mean_AD = mean(sub$EpilepsyScore[sub$Diagnosis == "AD"], na.rm = TRUE),
      Median_Control = median(sub$EpilepsyScore[sub$Diagnosis == "Control"], na.rm = TRUE),
      Median_AD = median(sub$EpilepsyScore[sub$Diagnosis == "AD"], na.rm = TRUE),
      Wilcoxon_P = wt$p.value,
      Welch_t_P = tt$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}))

write.csv(
  region_sensitivity,
  file = file.path(dir_tables, "GSE5281_signature_size_region_wise_sensitivity.csv"),
  row.names = FALSE
)

# =========================================================
# Part B. MTG AD subgroup concordance across signature sizes
# =========================================================

mtg_assign_list <- lapply(names(robust_list), function(nm) {
  df <- robust_list[[nm]]$data
  mtg_ad <- subset(df, BrainRegion == "MTG" & Diagnosis == "AD")
  
  median_cut <- median(mtg_ad$EpilepsyScore, na.rm = TRUE)
  mtg_ad$HighLow <- ifelse(mtg_ad$EpilepsyScore >= median_cut, "High", "Low")
  
  out <- mtg_ad[, c("Sample", "HighLow")]
  colnames(out)[2] <- nm
  out
})

mtg_assign_df <- Reduce(
  function(x, y) merge(x, y, by = "Sample", all = TRUE),
  mtg_assign_list
)

write.csv(
  mtg_assign_df,
  file = file.path(dir_tables, "GSE5281_MTG_AD_subgroup_assignment_top50_top100_top200.csv"),
  row.names = FALSE
)

pairwise_concordance <- function(df, a, b) {
  sub <- df[, c("Sample", a, b)]
  colnames(sub) <- c("Sample", "A", "B")
  
  agreement <- mean(sub$A == sub$B)
  
  high_A <- sub$Sample[sub$A == "High"]
  high_B <- sub$Sample[sub$B == "High"]
  
  low_A <- sub$Sample[sub$A == "Low"]
  low_B <- sub$Sample[sub$B == "Low"]
  
  data.frame(
    Signature_1 = a,
    Signature_2 = b,
    Agreement_rate = agreement,
    High_group_Jaccard_index = length(intersect(high_A, high_B)) / length(union(high_A, high_B)),
    Low_group_Jaccard_index = length(intersect(low_A, low_B)) / length(union(low_A, low_B)),
    stringsAsFactors = FALSE
  )
}

mtg_concordance <- do.call(rbind, list(
  pairwise_concordance(mtg_assign_df, "Top50", "Top100"),
  pairwise_concordance(mtg_assign_df, "Top50", "Top200"),
  pairwise_concordance(mtg_assign_df, "Top100", "Top200")
))

write.csv(
  mtg_concordance,
  file = file.path(dir_tables, "GSE5281_MTG_AD_subgroup_assignment_concordance.csv"),
  row.names = FALSE
)

# =========================================================
# Part C. Centered singscore robustness
# =========================================================

top100_up_file <- file.path(dir_signature, "Epilepsy_Up_top100.txt")
top100_down_file <- file.path(dir_signature, "Epilepsy_Down_top100.txt")

up_sig <- read_signature(top100_up_file)
down_sig <- read_signature(top100_down_file)

up_use <- intersect(up_sig, rownames(expr_gene))
down_use <- intersect(down_sig, rownames(expr_gene))

rank_data <- singscore::rankGenes(expr_gene)

ss_center <- singscore::simpleScore(
  rankData = rank_data,
  upSet = up_use,
  downSet = down_use,
  centerScore = TRUE,
  knownDirection = TRUE
)

if (!"TotalScore" %in% colnames(ss_center)) {
  stop("singscore output does not contain TotalScore.")
}

score_sing <- data.frame(
  Sample = rownames(ss_center),
  EpilepsyScore_singscore = as.numeric(ss_center[, "TotalScore"]),
  stringsAsFactors = FALSE
)

compare_df <- merge(
  score_df[, c("Sample", "EpilepsyScore", "Diagnosis", "BrainRegion")],
  score_sing,
  by = "Sample"
)

colnames(compare_df)[colnames(compare_df) == "EpilepsyScore"] <- "EpilepsyScore_ssGSEA"

write.csv(
  compare_df,
  file = file.path(dir_tables, "GSE5281_ssGSEA_vs_singscore_scores.csv"),
  row.names = FALSE
)

# -------------------------
# Correlation between scoring methods
# -------------------------
cor_res <- cor.test(
  compare_df$EpilepsyScore_ssGSEA,
  compare_df$EpilepsyScore_singscore,
  method = "spearman",
  exact = FALSE
)

cor_table <- data.frame(
  Comparison = "ssGSEA-based EpilepsyScore vs centered singscore-based projection score",
  Spearman_rho = unname(cor_res$estimate),
  P_value = cor_res$p.value,
  N = nrow(compare_df),
  stringsAsFactors = FALSE
)

write.csv(
  cor_table,
  file = file.path(dir_tables, "GSE5281_ssGSEA_vs_singscore_correlation.csv"),
  row.names = FALSE
)

# -------------------------
# Direction consistency table
# -------------------------
group_summary <- function(df, score_col, subset_name) {
  dd <- df[!is.na(df[[score_col]]) & !is.na(df$Diagnosis), ]
  
  if (length(unique(dd$Diagnosis)) < 2) {
    return(NULL)
  }
  
  ctrl_vals <- dd[[score_col]][dd$Diagnosis == "Control"]
  ad_vals <- dd[[score_col]][dd$Diagnosis == "AD"]
  
  wt <- wilcox.test(dd[[score_col]] ~ dd$Diagnosis)
  tt <- t.test(dd[[score_col]] ~ dd$Diagnosis)
  
  data.frame(
    Subset = subset_name,
    Score_method = score_col,
    N_total = nrow(dd),
    N_Control = length(ctrl_vals),
    N_AD = length(ad_vals),
    Mean_Control = mean(ctrl_vals, na.rm = TRUE),
    Mean_AD = mean(ad_vals, na.rm = TRUE),
    Median_Control = median(ctrl_vals, na.rm = TRUE),
    Median_AD = median(ad_vals, na.rm = TRUE),
    Direction = ifelse(mean(ad_vals, na.rm = TRUE) < mean(ctrl_vals, na.rm = TRUE),
                       "AD lower", "AD higher"),
    Wilcoxon_P = wt$p.value,
    Welch_t_P = tt$p.value,
    stringsAsFactors = FALSE
  )
}

summary_list <- list(
  group_summary(compare_df, "EpilepsyScore_ssGSEA", "Overall"),
  group_summary(compare_df, "EpilepsyScore_singscore", "Overall")
)

for (rg in c("MTG", "EC", "HIP", "SFG")) {
  sub <- subset(compare_df, BrainRegion == rg)
  summary_list[[length(summary_list) + 1]] <- group_summary(sub, "EpilepsyScore_ssGSEA", rg)
  summary_list[[length(summary_list) + 1]] <- group_summary(sub, "EpilepsyScore_singscore", rg)
}

direction_table <- do.call(rbind, summary_list)

write.csv(
  direction_table,
  file = file.path(dir_tables, "GSE5281_singscore_direction_consistency.csv"),
  row.names = FALSE
)

# -------------------------
# Scatter plot source and figure
# -------------------------
write.csv(
  compare_df,
  file = file.path(dir_figdata, "Supplementary_Figure_S6_source_data.csv"),
  row.names = FALSE
)

p_s6 <- ggplot(
  compare_df,
  aes(x = EpilepsyScore_ssGSEA, y = EpilepsyScore_singscore, color = Diagnosis)
) +
  geom_point(size = 2.2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = "black") +
  theme_bw(base_size = 12) +
  labs(
    x = "EpilepsyScore (ssGSEA)",
    y = "EpilepsyScore (centered singscore)",
    color = "Diagnosis"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(
  filename = file.path(dir_results, "Supplementary_Figure_S6_ssGSEA_vs_singscore.png"),
  plot = p_s6,
  width = 7.2,
  height = 5.6,
  dpi = 300
)

message("Signature-size sensitivity and singscore robustness analyses completed.")
