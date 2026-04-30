# HADACA3 — Benchmark Pipeline Modules

This document describes the modular pipeline used in the HADACA3 benchmark for multi-omics deconvolution. A complete pipeline is defined by selecting **one method per module**. The four modules are applied sequentially: preprocessing → feature selection → deconvolution → integration.

---

## Module 1 — Preprocessing

Applied independently to RNA-seq and DNA methylation data before feature selection and deconvolution.

| Method | Modality | Description |
|--------|----------|-------------|
| `ppID` | RNA + DNAm | No preprocessing (identity). Serves as baseline. |
| `Scale` | RNA + DNAm | Column-sum normalization. Each sample is divided by its total count, making all samples sum to 1. |
| `LogNorm` | RNA + DNAm | Log-normalization via Seurat's LogNormalize function (scale factor 10⁴ + log transform). For bulk RNA and DNAm, an exponential transformation is applied post hoc to recover linear scale. |

---

## Module 2 — Feature Selection

Selects a subset of genes (RNA) or CpG probes (DNAm) most informative for cell-type discrimination.

### RNA-seq

| Method | Description |
|--------|-------------|
| `fsID` | No selection. All ~20,000 genes are passed to deconvolution. Baseline. |
| `Toastbulknbfs` | Top 1,000 marker genes via TOAST applied to bulk RNA reference. |
| `Toastvst` | Top 1,000 marker genes via TOAST on VST-transformed bulk reference. |
| `SCcluster` | Differential expression markers from clustered scRNA-seq (Seurat Wilcoxon test, adj. p < 0.05). |
| `scpseudobulk` | Pairwise t-test markers from scRNA-seq; top 20 genes per cell-type pair across datasets. |

### DNA Methylation

| Method | Description |
|--------|-------------|
| `fsID` | No selection. All ~23,000 CpG probes are passed to deconvolution. Baseline. |
| `Toastpercent` | Top 80% probes ranked by cell-type specificity via TOAST. |
| `mostmethylated` | Probes above the 75th percentile of methylation per cell type (biologically informed). |
| `maxdiscriminant` | Maximally discriminant non-overlapping probes per cell type (up to 100 per type). |
| `splsda` | Sparse PLS-DA on logit-transformed reference; top 1,000 probes over 2 components. |

---

## Module 3 — Deconvolution

All methods solve a variant of: **Y ≈ X·p**, with **p ≥ 0** and **Σp = 1**.

| Method | Principle | Key idea |
|--------|-----------|----------|
| `lm` | Ordinary least squares | Minimizes ‖Y − Xp‖² without constraints; negatives set to zero post hoc. |
| `nnls` | Non-negative least squares | Enforces non-negativity by construction (Lawson-Hanson algorithm). |
| `nnlslargeref` | NNLS with reference truncation | Iteratively removes one cell type if it reduces reconstruction RMSE. |
| `epic` | Constrained least squares | Simplex-constrained optimization with internal TPM normalization. |
| `RLR` | Robust linear regression | Down-weights outlier features via iteratively reweighted least squares (IRLS / EpiDISH). |
| `RLRpoisson` | RLR with Poisson weights | Extends RLR by weighting features inversely to their mean expression (Poisson variance). |
| `RLRnnls` | Ensemble RLR + NNLS | Runs both independently; selects the estimate with lower reconstruction RMSE. |

> **Recommended:** RLR-based methods consistently achieve the best performance. `RLRpoisson` further improves RNA-seq results. `epic` is not recommended for DNAm data (TPM normalization incompatible with beta values).

---

## Module 4 — Multi-omics Integration


### Early Integration
Combines RNA and DNAm **at the feature level** before deconvolution.

| Method | Principle | Key idea | Main limitation |
|--------|-----------|----------|-----------------|
| `concatnoscale` | Raw concatenation | Stacks RNA and DNAm features without transformation. | Sensitive to scale differences between modalities. |
| `concatscale` | Normalized concatenation | Concatenation + sample-wise Gaussian CDF normalization. | Ignores cross-modal correlations. |
| `omicade4bulk` | Latent linear embedding | Joint low-dimensional representation via Multiple Co-Inertia Analysis (MCIA). | Sensitive to non-linear cross-omic relationships. |
| `Kernel` | Non-linear kernel embedding | Joint representation via kernel PCA over combined modality-specific kernels (STATIS-UMKL). | Sensitive to kernel choice; computationally expensive. |
| `OT` | Optimal transport | Cross-modal alignment via entropic OT in VAE latent space (uniPort). | High computational cost; depends on distance metric quality. |

### Late Integration
Performs deconvolution **independently per modality**, then combines predictions.

| Method | Principle | Key idea | Main limitation |
|--------|-----------|----------|-----------------|
| `onlyRNA` | RNA only | Uses RNA-based estimates only. | Ignores DNAm signal. |
| `onlyDNAm` | DNAm only | Uses DNAm-based estimates only. | Ignores RNA signal. |
| `limean` | Uniform averaging | Averages RNA and DNAm proportion estimates with equal weights. | Assumes equal reliability of both modalities. |
| `limeanRMSE` | Error-weighted aggregation | Weights modalities inversely proportional to their reconstruction RMSE. | Sensitive to reconstruction bias. |
| `tunedJ` | Rule-based selective averaging | Uses DNAm only for tumor subtypes (basal/classical); averages otherwise. | Requires prior biological knowledge. |




---

## Resources

⚠️ This paper is currently under anonymous review. All links to external resources (competition platform, code repository, datasets, and paper) will be made available upon completion of the review process.

- 🔗 **Competition platform:** [Codabench](https://www.codabench.org/competitions/4714/)
- 🔗 **Pipeline (Nextflow):** Github URL
- 🔗 **Datasets:** Zenodo URL + NCBI GEO URL
- 📄 **Paper:** *On the Promises and Limits of Multimodal Integration for Deconvolution: The HADACA3 Benchmark*

