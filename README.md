# COVID-19 Severity Interaction eQTL (ieQTL) Analysis Pipeline

## Overview

This pipeline identifies **interaction eQTLs (ieQTLs)**-genetic variants whose effect on gene expression is modulated by COVID-19 disease severity. It tests the model:

```
Expression ~ β₁·Genotype + β₂·Severity + β₃·(Genotype × Severity) + Covariates + ε
```

The interaction term (β₃) captures severity-dependent genetic regulation, revealing SNPs that have different regulatory effects in severe vs. non-severe COVID-19 patients.

## Study Design

| Parameter | Value |
|-----------|-------|
| Total samples | 430 (after QC) |
| Non-severe (Class = 0) | 312 |
| Severe (Class = 1) | 118 |
| Covariates | 69 (60 PEER + 5 genotype PCs + 4 binary) |
| Genes tested | 16,523 (autosomal, protein-coding) |
| Variants | ~8.2M (biallelic SNPs, autosomes 1–22) |
| Interaction term | Mean-centered severity (binary) |
| Genome build | GRCh38 (Gencode v25) |
| Data type | Whole blood RNA-seq + WGS |

## Pipeline Steps

### Step 1-Genotypes to 430 Samples

Subsets the original samples PLINK2 pgen to 430 QC-passing samples using a keep list derived from the covariate file.

```bash
plink2 --pfile <PGEN_PREFIX> \
       --keep keep_430_samples.txt \
       --make-pgen \
       --out TCD_UU_430_ieqtl
```

**Script:** `step1_subset_pgen.sh`

### Step 2-Build Mean-Centered Interaction File

interaction file (severity: 0=non-severe, 1=severe) to 430 samples and applies mean-centering.

Mean-centering ensures:
- `b_g` = average SNP effect across both groups (comparable to standard eQTL betas)
- `b_gi` = differential SNP effect between groups
- `p_gi` = ieQTL significance (unchanged by centering)

**Script:** `Prepare_interactionfile.R`  
**Output:** `interaction_430_centered.txt`, `interaction_430_lookup.csv`

### Step 3-Build Expression BED for TensorQTL

Converts the INT-transformed expression matrix (from the eQTL pipeline) into TensorQTL BED format:
- Merges with Gencode v25 gene coordinates (TSS-based)
- Removes mitochondrial genes (15 genes) and X/Y chromosome genes (537 genes)
- Retains 16,523 autosomal genes
- Sorts by chromosome and position
- Compresses with bgzip and indexes with tabix

**Script:** `prepare_bed_file.R`  
**Output:** `expression_430_tensorqtl.bed.gz` + `.tbi`

### Step 4-Run ieQTL Mapping

Runs TensorQTL `cis.map_nominal` with the interaction term. Includes a 4-way sample alignment validation (genotype × expression × covariates × interaction) before execution.

**Script:** `main_ieqtl_430.py`  
**SLURM:** `submit_ieqtl.sh`

#### TensorQTL Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `maf_threshold_interaction` | 0.05 | Ensures ≥12 minor allele carriers in smaller group (n=118) |
| `run_eigenmt` | True | Efficient multiple testing correction per gene |
| `write_top` | True | Outputs top association per gene |
| `write_stats` | True | Outputs full summary statistics |
| cis-window | ±1 Mb | Default TensorQTL window around TSS |

## Input Files

| File | Description | Dimensions |
|------|-------------|------------|
| `TCD_UU_430_ieqtl.{pgen,pvar,psam}` | 430-sample genotype files | ~8.2M variants × 430 samples |
| `expression_430_tensorqtl.bed.gz` | INT-normalised expression | 16,523 genes × 430 samples |
| `<COVARIATE_FILE_430>` | Covariates (PEER + PCs + binary) | 69 covariates × 430 samples |
| `interaction_430_centered.txt` | Mean-centered severity term | 430 samples × 1 term |

## Output Files

| File | Description |
|------|-------------|
| `<PREFIX>.cis_qtl_top_assoc.txt.gz` | Top SNP per gene with ieQTL statistics |
| `<PREFIX>.cis_qtl_pairs.<chr>.parquet` | Full nominal results per chromosome |

### Key Output Columns

| Column | Description |
|--------|-------------|
| `b_g` | Average SNP effect on expression (due to mean-centering) |
| `b_i` | Severity main effect on expression |
| `b_gi` | **Interaction effect**-change in SNP effect between severity groups |
| `pval_gi` | **ieQTL p-value**-significance of the interaction term |
| `pval_adj_bh` | BH-adjusted p-value (if available) |

## Expression Data Processing (Upstream)

The expression BED used here is the end product of the eQTL preprocessing pipeline:

1. Raw count merging across two cohorts (AA: 320, BB: 110 samples)
2. Gene ID cleaning (Ensembl version stripping, PAR gene deduplication)
3. GTEx-style dual filtering (TPM ≥ 0.1 AND counts ≥ 6 in ≥20% samples)
4. ComBat-seq batch correction (cohort-level systematic differences)
5. TMM normalisation (edgeR)
6. Inverse Normal Transformation (INT) per gene
7. 60 PEER factors computed with `PEER_setAdd_mean=TRUE`

## Scientific Design Decisions

### Why binary severity (not 3-group)?

The original clinical classification included mild, moderate, and severe groups (~138 / ~176 / ~119). We merged mild + moderate into non-severe because:

- **Power:** N=430 is in the lower range for ieQTL detection; splitting further reduces power in each stratum
- **Linearity assumption:** TensorQTL treats the interaction term as continuous-a 0/1/2 coding assumes linear dose-response, which may not hold
- **Comparability:** Wang et al. (2022) used binary severity with a comparable sample size (n=465; 359 severe, 106 non-severe)
- **Validation strategy:** 3-group gradient analysis on top hits serves as supplementary validation

### Why mean-center the interaction term?

With raw 0/1 coding, `b_g` represents the SNP effect only in the reference group (non-severe). Mean-centering makes `b_g` the average effect across groups, which is more interpretable and directly comparable to standard eQTL results. The ieQTL p-value (`pval_gi`) is invariant to centering.

### Why exclude severity from covariates?

Severity is the study variable of interest. Including it as a covariate would regress out the biological signal being tested by the interaction term.

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| [TensorQTL](https://github.com/broadinstitute/tensorqtl) | ≥1.0.7 | ieQTL mapping (GPU-accelerated) |
| [PLINK2](https://www.cog-genomics.org/plink/2.0/) | ≥2.0 | Genotype subsetting |
| [htslib](http://www.htslib.org/) | ≥1.9 | bgzip + tabix |
| R ≥ 4.0 |-| Data preparation |
| Python ≥ 3.8 |-| TensorQTL execution |
| PyTorch | ≥1.7 | GPU/CPU tensor operations |
| CUDA (optional) | ≥11.0 | GPU acceleration |

### R Packages

`data.table`

### Python Packages

`pandas`, `torch`, `tensorqtl`, `scipy`

## References

- Wang QS, Edahiro R, Namkoong H, et al. The whole blood transcriptional regulation landscape in 465 COVID-19 infected samples from Japan COVID-19 Task Force. *Nature Communications* 13, 4830 (2022). [DOI: 10.1038/s41467-022-32276-2](https://doi.org/10.1038/s41467-022-32276-2)
- Taylor-Weiner A, Aguet F, et al. Scaling computational genomics to millions of individuals with GPUs. *Genome Biology* 20, 228 (2019).
- Ongen H, Buil A, Brown AA, Dermitzakis ET, Delaneau O. Fast and efficient QTL mapper for thousands of molecular phenotypes. *Bioinformatics* 32, 1479–1485 (2016).

## License

[Specify your license here]

## Contact

[Your contact information]
