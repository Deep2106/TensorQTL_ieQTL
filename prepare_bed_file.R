# ============================================================
# Step 3: Build TensorQTL-format expression BED from
#         INT-transformed expression + gene coordinates
#         Remove MT genes, filter to autosomes, no chr prefix
# ============================================================

library(data.table)

WORK_DIR   <- "/path/to/ieQTL"
EXPR_FILE  <- "expr_int_final_430.txt"              # expr_int_final_430.txt
COORD_FILE <- "gene_positions_matrixqtl.txt"                 # gene_positions_matrixqtl.txt
PSAM_FILE  <- file.path(WORK_DIR, "TCD_UU_430_ieqtl.psam")
OUT_BED    <- file.path(WORK_DIR, "expression_430_tensorqtl.bed")

# --- Load expression ---
expr <- fread(EXPR_FILE)
gene_ids <- expr$ids
expr_mat <- expr[, -1]
cat("Expression loaded:", nrow(expr), "genes x", ncol(expr_mat), "samples\n")

# --- Load gene coordinates ---
coords <- fread(COORD_FILE)
cat("Coordinates loaded:", nrow(coords), "genes\n")

# --- Filter coordinates to autosomes only (chr 1-22) ---
coords <- coords[CHR %in% as.character(1:22), ]
cat("After autosome filter:", nrow(coords), "genes\n")

# --- Identify genes to remove ---
# MT genes: in expression but not in any coordinates
mt_genes <- gene_ids[!gene_ids %in% fread(COORD_FILE)$gene_id]
cat("\nMT genes (no coordinates):", length(mt_genes), "\n")
cat(mt_genes, sep = "\n")

# X/Y genes: in coordinates but removed by autosome filter
xy_genes <- fread(COORD_FILE)[CHR %in% c("X", "Y"), gene_id]
cat("\nX/Y genes removed:", length(xy_genes), "\n")

# --- Filter expression to autosomal genes with coordinates ---
keep_genes <- gene_ids %in% coords$gene_id
expr_filt <- expr_mat[keep_genes, ]
gene_ids_filt <- gene_ids[keep_genes]
cat("\nAfter filtering (MT + X/Y removed):", length(gene_ids_filt), "genes\n")
stopifnot(length(gene_ids_filt) == nrow(coords))

# --- Align coordinates to expression gene order ---
coord_order <- match(gene_ids_filt, coords$gene_id)
coords_aligned <- coords[coord_order, ]
stopifnot(all(coords_aligned$gene_id == gene_ids_filt))

# --- Build BED data frame ---
# TensorQTL format: #chr  start  end  gene_id  sample1  sample2 ...
# No chr prefix - plain integers (1, 2, ..., 22)
bed <- data.frame(
  `#chr`    = coords_aligned$CHR,
  start     = coords_aligned$TSS,
  end       = coords_aligned$end,
  gene_id   = coords_aligned$gene_id,
  check.names = FALSE
)
bed <- cbind(bed, expr_filt)

# --- Sort by chromosome (numeric) then by start position ---
chr_num <- as.integer(coords_aligned$CHR)
sort_idx <- order(chr_num, coords_aligned$TSS)
bed <- bed[sort_idx, ]

cat("\n=== Final BED summary ===\n")
cat("Dimensions:", nrow(bed), "genes x", ncol(bed) - 4, "samples\n")
cat("Chromosomes:", sort(unique(as.integer(bed$`#chr`))), "\n")
cat("\nFirst 3 rows (first 6 cols):\n")
print(bed[1:3, 1:6])
cat("\nLast 3 rows (first 6 cols):\n")
print(bed[(nrow(bed)-2):nrow(bed), 1:6])

# --- Write BED ---
fwrite(bed, OUT_BED, sep = "\t", quote = FALSE)
cat("\nBED written to:", OUT_BED, "\n")

# --- Cross-check sample IDs with pgen psam ---
psam <- fread(PSAM_FILE)
psam_ids <- as.character(psam$IID)
bed_sample_ids <- colnames(bed)[5:ncol(bed)]

cat("\n=== Sample ID cross-check (BED vs PSAM) ===\n")
cat("BED samples:", length(bed_sample_ids), "\n")
cat("PSAM samples:", length(psam_ids), "\n")
cat("Overlap:", sum(bed_sample_ids %in% psam_ids), "\n")
cat("In BED but not PSAM:", sum(!bed_sample_ids %in% psam_ids), "\n")
cat("In PSAM but not BED:", sum(!psam_ids %in% bed_sample_ids), "\n")

if (all(bed_sample_ids %in% psam_ids) & all(psam_ids %in% bed_sample_ids)) {
  cat("PASS: All 430 sample IDs match perfectly\n")
} else {
  cat("FAIL: Sample ID mismatch detected - investigate before proceeding\n")
}

# --- Cross-check sample IDs with interaction file ---
INT_FILE <- file.path(WORK_DIR, "interaction_430_centered.txt")
int_df <- fread(INT_FILE)
int_ids <- as.character(int_df$`#IDS`)

cat("\n=== Sample ID cross-check (BED vs Interaction) ===\n")
cat("Overlap:", sum(bed_sample_ids %in% int_ids), "\n")
cat("In BED but not Interaction:", sum(!bed_sample_ids %in% int_ids), "\n")
cat("In Interaction but not BED:", sum(!int_ids %in% bed_sample_ids), "\n")

if (all(bed_sample_ids %in% int_ids) & all(int_ids %in% bed_sample_ids)) {
  cat("PASS: All 430 sample IDs match perfectly\n")
} else {
  cat("FAIL: Sample ID mismatch detected - investigate before proceeding\n")
}

# --- Cross-check sample IDs with covariate file ---
COV_FILE <- "Covariates_TCD_UU_Peer60_latest.tsv"
cov_ids <- colnames(fread(COV_FILE, nrows = 0))[-1]

cat("\n=== Sample ID cross-check (BED vs Covariates) ===\n")
cat("Overlap:", sum(bed_sample_ids %in% cov_ids), "\n")

if (all(bed_sample_ids %in% cov_ids) & all(cov_ids %in% bed_sample_ids)) {
  cat("PASS: All 430 sample IDs match perfectly\n")
} else {
  cat("FAIL: Sample ID mismatch detected - investigate before proceeding\n")
}
