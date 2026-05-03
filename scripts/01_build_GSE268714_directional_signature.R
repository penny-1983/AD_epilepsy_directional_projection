# =========================================================
# 01_build_GSE268714_directional_signature.R
# Purpose:
# Rebuild direction-aware epilepsy-derived SOZ-up and SOZ-down signatures
# from GSE268714 bulk RNA-seq count data.
# =========================================================

source("scripts/00_setup.R")

# -------------------------
# Required packages
# -------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

# -------------------------
# Input file
# -------------------------
count_file <- file.path(dir_raw, "GSE268714", "GSE268714_EPILEPSY.BULK.RNAseq.Count.txt")

if (!file.exists(count_file)) {
  stop("Count file not found. Please download GSE268714 count data from GEO and place it in: ",
       file.path(dir_raw, "GSE268714"))
}

# -------------------------
# Read count matrix
# -------------------------
count_df <- read.delim(
  count_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_col <- colnames(count_df)[1]

count_mat <- as.data.frame(count_df)
rownames(count_mat) <- count_mat[[gene_col]]
count_mat <- count_mat[, -1, drop = FALSE]
count_mat <- as.matrix(count_mat)
mode(count_mat) <- "numeric"

# -------------------------
# Build sample group annotation
# Cold  -> NIZ
# Onset -> SOZ
# Spread -> PZ
# -------------------------
sample_names <- colnames(count_mat)

group_full <- ifelse(
  grepl("Cold", sample_names, ignore.case = TRUE), "NIZ",
  ifelse(
    grepl("Onset", sample_names, ignore.case = TRUE), "SOZ",
    ifelse(grepl("Spread", sample_names, ignore.case = TRUE), "PZ", NA)
  )
)

group_info <- data.frame(
  Sample = sample_names,
  Group = group_full,
  stringsAsFactors = FALSE
)

write.csv(
  group_info,
  file = file.path(dir_tables, "GSE268714_sample_group_mapping.csv"),
  row.names = FALSE
)

# -------------------------
# Keep NIZ and SOZ samples for primary contrast
# -------------------------
keep_samples <- group_info$Group %in% c("NIZ", "SOZ")
count_sub <- count_mat[, keep_samples, drop = FALSE]
group <- factor(group_info$Group[keep_samples], levels = c("NIZ", "SOZ"))

# -------------------------
# Differential expression analysis
# -------------------------
dge <- DGEList(counts = count_sub, group = group)
keep_gene <- filterByExpr(dge, group = group)
dge <- dge[keep_gene, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ group)

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

deg_res <- topTable(fit, coef = "groupSOZ", number = Inf, sort.by = "none")
deg_res$Gene <- rownames(deg_res)

deg_res <- deg_res[, c("Gene", setdiff(colnames(deg_res), "Gene"))]

write.csv(
  deg_res,
  file = file.path(dir_tables, "GSE268714_SOZ_vs_NIZ_limma_results.csv"),
  row.names = FALSE
)

# -------------------------
# Direction-aware signatures
# -------------------------
up_res <- deg_res[deg_res$logFC > 0, ]
down_res <- deg_res[deg_res$logFC < 0, ]

up_res <- up_res[order(-up_res$logFC), ]
down_res <- down_res[order(down_res$logFC), ]

signature_sizes <- c(50, 100, 200)

for (n in signature_sizes) {
  up_genes <- head(up_res$Gene, n)
  down_genes <- head(down_res$Gene, n)
  
  writeLines(
    up_genes,
    con = file.path(dir_signature, paste0("Epilepsy_Up_top", n, ".txt"))
  )
  
  writeLines(
    down_genes,
    con = file.path(dir_signature, paste0("Epilepsy_Down_top", n, ".txt"))
  )
}

# -------------------------
# Export top100 combined signature table
# -------------------------
top100_up <- head(up_res, 100)
top100_down <- head(down_res, 100)

top100_up$Direction <- "SOZ_up"
top100_down$Direction <- "SOZ_down"

sig100 <- rbind(top100_up, top100_down)

write.csv(
  sig100,
  file = file.path(dir_tables, "GSE268714_directional_signature_top100.csv"),
  row.names = FALSE
)

# -------------------------
# Save normalization summary
# -------------------------
dge_summary <- as.data.frame(dge$samples)
dge_summary$Sample <- rownames(dge_summary)

write.csv(
  dge_summary,
  file = file.path(dir_tables, "GSE268714_TMM_normalization_summary.csv"),
  row.names = FALSE
)

message("GSE268714 directional signature construction completed.")
