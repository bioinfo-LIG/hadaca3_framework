# 00_paper_reproducibility

This folder contains all the code and resources needed to reproduce the results, figures, and analyses presented in the paper:

> *On the Promises and Limits of Multimodal Integration for Deconvolution: The HADACA3 Benchmark* 

> ⚠️ **This paper is currently under anonymous review. All links to external resources will be made available upon completion of the review process.**

---


## Folder Structure

```
00_paper_reproducibility/
├── 1_data_simulations/       # Scripts to generate the in silico benchmark datasets
├── 2_scripts_best_pipelines/ # Scripts reproducing the best-performing pipelines
└── 3_figures_code/           # Code to reproduce all figures in the paper
```

---

## 1 — `1_data_simulations/`


The following simulation scenarios are covered:

| Dataset | Noise model (RNA) | Noise model (DNAm) | Notes |
|---------|-------------------|---------------------|-------|
| `SBN5` | Pseudo-bulk | Pseudo-bulk | Aggregated single-cell profiles |
| `SDN5` | χ² (heteroscedastic) | Gaussian | 5 cell types, matching reference |
| `SDN4` | χ² (heteroscedastic) | Gaussian | 4 cell types (basal-like removed) |
| `SDN6` | χ² (heteroscedastic) | Gaussian | 6 cell types (extra synthetic type) |
| `SDE5` | EM (low-rank factor) | EM (low-rank factor) | Structured gene co-expression |
| `SDEL` | EM (low-rank factor) | EM (low-rank factor) | Rare cell type scenario |
| `SDC5` | Copula + Negative Binomial | Copula + Beta | Structured cross-feature dependencies |



---

## 2 — `2_scripts_best_pipelines/`

Contains the scripts reproducing the **best-performing pipeline combinations** identified in the benchmark, for both early and late integration strategies. 



---

## 3 — `3_figures_code/`

Contains the code to reproduce main **figures** presented in the main paper and supplementary material.

