# uniPort early integration (block version 1.0, 2026-11-13).
#
# Diagonal integration with uniPort (Cao et al., Nat Commun 2022): one variational encoder per omic (RNA, DNAm),
# the two latent spaces aligned by unbalanced optimal transport (uniPort mode 'd'). Each sample - a mixture or a
# reference cell type - has an RNA and a DNAm encoding: its latent profile is the mean of the two (the samples are
# paired across the omics). As for MCIA and MKL, the latent coordinates are shifted to be non-negative (minus the
# global minimum), and the latent profiles of the mixtures and of the reference cell types form the integrated unit
# that is passed to the deconvolution.
#
# This block replaces function_blocks/early_integration/OT.py (block "OT"), which ran uniPort in mode 'v' and returned
# adata.X, i.e. the normalised RNA data: its latent space was never read, the DNAm did not enter the prediction, and
# mode 'v' has no optimal-transport term. The block "OT" is kept, unchanged and annotated, for the record only.
#
# Data preparation: as in the block "OT" (per omic: normalize_total, log1p, division by the global maximum; samples
# with an all-zero profile in either omic removed). Settings: 500 iterations, batch size min(50, number of samples),
# 16 latent dimensions (uniPort default), random seed 124 (uniPort default, set explicitly here), 2 torch threads.
#
# Reproducibility: with the same seed, the result depends on the number of torch threads (and may depend on the
# platform and on the versions of torch / BLAS): on three configurations, 1 thread instead of 2 changed the
# integration by up to 0.45 and the aggregated scores by -0.11 to +0.03 (one pipeline from 0.53 to 0). The number of
# threads is therefore fixed here (THREADS), not left to the environment. Results are reproducible for a given
# platform and software versions; across platforms, expect run-to-run variation of this order.
#
# Software used for the first run: uniport 1.3, torch 2.14.1 (CPU), scanpy (see the MANIFEST of the run).

import os
import tempfile
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import uniport as up

BLOCK_VERSION = "1.0"
SEED = 124
ITERATIONS = 500
THREADS = 2  # the result depends on it (see above)


def program_block_EI(rna_unit, met_unit, path_og_dataset=""):
    torch.set_num_threads(THREADS)

    # ---- data preparation (as in the block "OT") ----
    common_genes = rna_unit["mix"].index.intersection(rna_unit["ref"].index)
    rna_mix, rna_ref = rna_unit["mix"].loc[common_genes], rna_unit["ref"].loc[common_genes]
    common_met = met_unit["mix"].index.intersection(met_unit["ref"].index)
    met_mix, met_ref = met_unit["mix"].loc[common_met], met_unit["ref"].loc[common_met]
    if list(rna_mix.columns) != list(met_mix.columns):
        met_mix = met_mix[rna_mix.columns]  # same mixture samples, same order
    met_ref = met_ref[rna_ref.columns]      # same cell types, same order

    n_sample, n_ref = rna_mix.shape[1], rna_ref.shape[1]
    types = ["sample"] * n_sample + ["ref"] * n_ref

    def make_adata(mix, ref, source, domain):
        a = sc.AnnData(pd.concat([mix, ref], axis=1).T.astype(np.float32))
        a.obs["source"] = source
        a.obs["type"] = types
        a.obs["domain_id"] = domain
        sc.pp.normalize_total(a)
        sc.pp.log1p(a)
        a.X = a.X / a.X.max()
        a.X = np.nan_to_num(a.X)
        return a

    ad_rna = make_adata(rna_mix, rna_ref, "RNA", 0)
    ad_met = make_adata(met_mix, met_ref, "MET", 1)
    keep = (np.asarray(ad_rna.X.sum(axis=1)).ravel() > 0) & (np.asarray(ad_met.X.sum(axis=1)).ravel() > 0)
    ad_rna, ad_met = ad_rna[keep].copy(), ad_met[keep].copy()
    if ad_rna.n_obs == 0:
        raise ValueError("no valid sample after filtering")

    # ---- integration: uniPort mode 'd' (one encoder per omic, latent spaces aligned by unbalanced OT) ----
    outdir = tempfile.mkdtemp(prefix="uniport_")  # own checkpoint folder (runs in parallel)
    up.Run(adatas=[ad_rna, ad_met], mode="d", iteration=ITERATIONS, batch_size=min(50, ad_rna.n_obs),
           outdir=outdir, seed=SEED)
    latent = (np.asarray(ad_rna.obsm["latent"]) + np.asarray(ad_met.obsm["latent"])) / 2
    latent = latent - latent.min()

    z = pd.DataFrame(latent.T, index=["L%d" % (i + 1) for i in range(latent.shape[1])], columns=ad_rna.obs_names)
    is_ref = (ad_rna.obs["type"] == "ref").values
    return {"mix": z.loc[:, ~is_ref], "ref": z.loc[:, is_ref], "ref_scRNA": rna_unit["ref_scRNA"]}
