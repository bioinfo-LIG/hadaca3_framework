# data_creation

This folder contains all R scripts used to generate the **reference profiles** and **benchmark datasets** described in the HADACA3 study. It is part of the `00_paper_reproducibility` directory.

> ⚠️ **This paper is currently under anonymous review. Citation and repository links will be made available upon completion of the review process.**


---

## 01.1 — Reference Profiles

Scripts in this folder generate the **cell-type reference profiles** provided to participants for supervised deconvolution. Three types of reference profiles are produced:

| Reference type | Modality | Cell types |
|----------------|----------|------------|
| Pure-bulk RNA | RNA-seq | Classical, Basal-like, Immune, Fibroblasts, Endothelial |
| Pure-bulk DNAm | DNA methylation | Classical, Basal-like, Immune, Fibroblasts, Endothelial |
| scRNA-seq | Single-cell RNA | Classical, Basal-like, Immune, Fibroblasts, Endothelial |

For bulk references, replicates from each cell type are aggregated into a single meta-reference profile by averaging all available cell-type-specific profiles.

> **All source data files are publicly available.** They can be downloaded from the repositories referenced in Tables 5 and 6 of the paper's Supplementary Material (Section B.1.1). Each source dataset is identified by its GEO accession number or public repository URL as cited in the paper.

---

## 01.2 — Pure Profiles for Simulation

Scripts in this folder generate the **pure cell-type profiles** used as input to the simulation models. These profiles are distinct from the reference profiles provided to participants, ensuring that the benchmark datasets do not overlap with the deconvolution references (no circularity).

> **All source data files are publicly available** and can be downloaded from the repositories referenced in Tables 5 and 6 of the paper's Supplementary Material (Section B.1.1).

---

## 01.2 — Simulation Scripts

The following datasets are generated:

| Dataset | Script | Noise RNA | Noise DNAm | Cell types | Notes |
|---------|--------|-----------|------------|------------|-------|
| `SBN5` | `sim_pseudobulk.R` | Pseudo-bulk | Pseudo-bulk | 5 | Aggregated single-cell profiles |
| `SDN5` | `sim_noDep_5CT.R` | χ² (heteroscedastic) | Gaussian | 5 | Matching reference |
| `SDN4` | `sim_noDep_4CT.R` | χ² (heteroscedastic) | Gaussian | 4 | Basal-like type removed |
| `SDN6` | `sim_noDep_6CT.R` | χ² (heteroscedastic) | Gaussian | 6 | Extra synthetic cell type |
| `SDE5` | `sim_EMFA.R` | EM (low-rank factor) | EM (low-rank factor) | 5 | Structured gene co-expression |
| `SDEL` | `sim_EMFA_rare.R` | EM (low-rank factor) | EM (low-rank factor) | 5 | Rare cell type scenario (`α_rare`) |
| `SDC5` | `sim_copula.R` | Copula + Negative Binomial | Copula + Beta | 5 | Structured cross-feature dependencies |

Each script produces two splits of 30 samples each (train and test), stored in HDF5 format.

---

## Requirements

- R ≥ 4.0
- R packages: `rhdf5`, `MASS`, `copula`, `Matrix`


---

## Notes

- All simulation parameters are documented within each script.
- The datasets are also available on Zenodo (DOI available after review).
