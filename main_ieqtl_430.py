# ============================================================
# ieQTL Analysis using TensorQTL
# Model: Expression ~ β1·SNP + β2·Severity + β3·(SNP×Severity) + Covariates
# Interaction term: mean-centered severity (0/1 → centered)
# 430 samples, 16,523 autosomal genes, ~8.2M variants
# ============================================================

import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis
import os

# --- Device setup ---
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")
print(f"tensorqtl: {tensorqtl.__version__}")

# --- Define paths (EDIT THESE) ---
WORK_DIR          = "/path/to/ieQTL/LATEST_IEQTL_RUN_March2026"
PGEN_PREFIX       = os.path.join(WORK_DIR, "TCD_UU_430_ieqtl")
EXPRESSION_BED    = os.path.join(WORK_DIR, "expression_430_tensorqtl.bed.gz")
COVARIATES_FILE   = "Covariates_TCD_UU_Peer60_latest.tsv"
INTERACTION_FILE  = os.path.join(WORK_DIR, "interaction_430_centered.txt")
OUTPUT_DIR        = os.path.join(WORK_DIR, "ieqtl_results")
PREFIX            = "Covid_430_ieqtl_severity_peer60"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Step 1: Load all inputs
# ============================================================

# --- Phenotypes (expression BED) ---
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(EXPRESSION_BED)
print(f"\nPhenotypes: {phenotype_df.shape[0]} genes x {phenotype_df.shape[1]} samples")

# --- Covariates ---
covariates_df = pd.read_csv(COVARIATES_FILE, sep='\t', index_col=0).T
print(f"Covariates: {covariates_df.shape[0]} samples x {covariates_df.shape[1]} covariates")

# --- Interaction term (mean-centered severity) ---
interaction_df = pd.read_csv(INTERACTION_FILE, sep='\t', index_col=0)
print(f"Interaction: {interaction_df.shape[0]} samples x {interaction_df.shape[1]} term(s)")
print(f"  Mean of interaction term: {interaction_df['Class'].mean():.10f} (should be ~0)")
print(f"  Unique values: {sorted(interaction_df['Class'].unique())}")

# --- Genotypes (pgen) ---
pgr = pgen.PgenReader(PGEN_PREFIX)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df
print(f"Genotypes: {genotype_df.shape[0]} variants x {genotype_df.shape[1]} samples")

# ============================================================
# Step 2: Validate sample alignment across all four matrices
# ============================================================

pheno_samples = set(phenotype_df.columns)
cov_samples   = set(covariates_df.index)
int_samples   = set(interaction_df.index)
geno_samples  = set(genotype_df.columns)

print(f"\n=== Sample Alignment Check ===")
print(f"Phenotype samples:   {len(pheno_samples)}")
print(f"Covariate samples:   {len(cov_samples)}")
print(f"Interaction samples: {len(int_samples)}")
print(f"Genotype samples:    {len(geno_samples)}")

# All pairwise intersections should be 430
intersection = pheno_samples & cov_samples & int_samples & geno_samples
print(f"4-way intersection:  {len(intersection)}")

assert len(intersection) == 430, \
    f"FATAL: Expected 430 samples in intersection, got {len(intersection)}"

# Check for mismatches
for name, s in [("Phenotype", pheno_samples), ("Covariate", cov_samples),
                ("Interaction", int_samples), ("Genotype", geno_samples)]:
    diff = s - intersection
    if diff:
        print(f"  WARNING: {name} has {len(diff)} extra samples: {diff}")

print("PASS: All 430 samples aligned across all four matrices")

# ============================================================
# Step 3: Run ieQTL mapping (cis nominal with interaction)
# ============================================================

print(f"\n=== Starting ieQTL Mapping ===")
print(f"  Genes:       {phenotype_df.shape[0]}")
print(f"  Variants:    {genotype_df.shape[0]}")
print(f"  Covariates:  {covariates_df.shape[1]}")
print(f"  Interaction: severity (mean-centered)")
print(f"  MAF threshold (interaction): 0.05")
print(f"  Output dir:  {OUTPUT_DIR}")
print(f"  Prefix:      {PREFIX}")

cis.map_nominal(
    genotype_df,
    variant_df,
    phenotype_df,
    phenotype_pos_df,
    PREFIX,
    covariates_df=covariates_df,
    interaction_df=interaction_df,
    maf_threshold_interaction=0.05,
    run_eigenmt=True,
    output_dir=OUTPUT_DIR,
    write_top=True,
    write_stats=True
)

print(f"\n=== ieQTL Mapping Complete ===")
print(f"Results in: {OUTPUT_DIR}")

# ============================================================
# Step 4: Quick summary of results
# ============================================================

# Load top associations (write_top=True produces this file)
top_file = os.path.join(OUTPUT_DIR, f"{PREFIX}.cis_qtl_top_assoc.txt.gz")
if os.path.exists(top_file):
    top_df = pd.read_csv(top_file, sep='\t')
    print(f"\n=== Top Associations Summary ===")
    print(f"Total genes tested: {len(top_df)}")

    # ieQTL significance (interaction p-value)
    n_sig_nom = (top_df['pval_gi'] < 5e-8).sum()
    n_sig_fdr = (top_df['pval_adj_bh'] < 0.05).sum() if 'pval_adj_bh' in top_df.columns else 'N/A'
    print(f"Significant ieQTLs (p_gi < 5e-8):  {n_sig_nom}")
    print(f"Significant ieQTLs (FDR < 0.05):   {n_sig_fdr}")

    # Show top 10 by interaction p-value
    print(f"\nTop 10 ieQTLs by interaction p-value:")
    top10 = top_df.nsmallest(10, 'pval_gi')[
        ['phenotype_id', 'variant_id', 'pval_gi', 'b_g', 'b_i', 'b_gi']
    ]
    print(top10.to_string(index=False))
else:
    print(f"Top associations file not found at {top_file}")
