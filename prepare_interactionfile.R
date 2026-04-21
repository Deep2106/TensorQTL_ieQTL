# ============================================================
# Step 2: Prepare interaction file and mean-center the interaction term
# ============================================================

library(data.table)

WORK_DIR    <- "/path/to/ieQTL"
INT_FILE    <- "interaction_2grp.txt"
COV_FILE    <- "Covariates_TCD_UU_Peer60_latest.tsv"
OUT_FILE    <- file.path(WORK_DIR, "interaction_430_centered.txt")

# --- Load interaction file ---
int_df <- fread(INT_FILE)
colnames(int_df) <- c("IDS", "Class")
int_df$IDS <- as.character(int_df$IDS)
cat("Interaction file loaded:", nrow(int_df), "samples\n")

# --- Get 430 sample IDs from covariate header ---
cov_header <- colnames(fread(COV_FILE, nrows = 0))
keep_ids   <- cov_header[-1]
cat("Samples from covariate file:", length(keep_ids), "\n")
stopifnot(length(keep_ids) == 430)

# --- Subset ---
int_430 <- int_df[int_df$IDS %in% keep_ids, ]
cat("After subsetting:", nrow(int_430), "samples\n")
stopifnot(nrow(int_430) == 430)

# --- Report raw split ---
cat("\n=== RAW severity split (430 samples) ===\n")
print(table(int_430$Class))

# --- Mean-center the interaction term ---
# Raw:      non-severe = 0, severe = 1
# Centered: non-severe = -mean, severe = 1-mean
# This ensures b_g = average genetic effect across groups
raw_mean <- mean(int_430$Class)
int_430$Class_centered <- int_430$Class - raw_mean

cat("\n=== Mean-centering summary ===\n")
cat("Raw mean (proportion severe):", round(raw_mean, 4), "\n")
cat("Non-severe coded as:", round(-raw_mean, 4), "\n")
cat("Severe coded as:    ", round(1 - raw_mean, 4), "\n")
cat("Centered mean (should be ~0):", round(mean(int_430$Class_centered), 10), "\n")

# --- Sanity check: verify centering ---
stopifnot(abs(mean(int_430$Class_centered)) < 1e-10)

# --- Write output ---
# TensorQTL expects: column 1 = sample IDs (index), column 2+ = interaction term(s)
# We write the centered version only
out_df <- data.frame(
  IDS = int_430$IDS,
  Class = int_430$Class_centered
)
colnames(out_df)[1] <- "#IDS"
fwrite(out_df, OUT_FILE, sep = "\t", quote = FALSE)
cat("\nWritten to:", OUT_FILE, "\n")

# --- Also save a lookup of raw → centered for reference ---
lookup <- data.frame(
  IDS = int_430$IDS,
  Class_raw = int_430$Class,
  Class_centered = int_430$Class_centered
)
fwrite(lookup, file.path(WORK_DIR, "interaction_430_lookup.csv"))
cat("Lookup table saved to: interaction_430_lookup.csv\n")
