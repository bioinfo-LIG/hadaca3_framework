
import uniport as up
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix


def program_block_EI(rna_unit,met_unit,path_og_dataset=""):

    # Align by shared genes
    common_genes = rna_unit["mix"].index.intersection(rna_unit["ref"].index)
    rna_unit["mix"] = rna_unit["mix"].loc[common_genes]
    rna_unit["ref"] = rna_unit["ref"].loc[common_genes]

    # Rebuild combined dataset
    data_rna_DT = pd.concat([rna_unit["mix"], rna_unit["ref"]], axis=1)
    data_rna_D = rna_unit["mix"]

    # Do the same for MET if needed
    common_genes_met = met_unit["mix"].index.intersection(met_unit["ref"].index)
    met_unit["mix"] = met_unit["mix"].loc[common_genes_met]
    met_unit["ref"] = met_unit["ref"].loc[common_genes_met]
    data_met_DT = pd.concat([met_unit["mix"], met_unit["ref"]], axis=1)

    # Make anndata
    n_sample = data_rna_D.shape[1]
    n_ref = data_rna_DT.shape[1] - n_sample
    adata_rna_DT = sc.AnnData(data_rna_DT.T, dtype=np.float32)
    adata_met_DT = sc.AnnData(data_met_DT.T, dtype=np.float32)
    adata_rna_DT.obs['source'] = 'RNA'
    adata_rna_DT.obs['type'] = ['sample']*n_sample + ['ref']*n_ref
    adata_met_DT.obs['source'] = 'MET'
    adata_met_DT.obs['type'] = ['sample']*n_sample + ['ref']*n_ref


    # # Select 2k HVG in each
    sc.pp.normalize_total(adata_rna_DT)
    sc.pp.log1p(adata_rna_DT)

    sc.pp.normalize_total(adata_met_DT)
    sc.pp.log1p(adata_met_DT)
    adata_rna_DT.X = adata_rna_DT.X/adata_rna_DT.X.max()
    adata_met_DT.X = adata_met_DT.X/adata_met_DT.X.max()


    # --- Safety checks ---
    adata_rna_DT.X = np.nan_to_num(adata_rna_DT.X)
    adata_met_DT.X = np.nan_to_num(adata_met_DT.X)

    adata_rna_DT = adata_rna_DT[adata_rna_DT.X.sum(axis=1) > 0]
    adata_met_DT = adata_met_DT[adata_met_DT.X.sum(axis=1) > 0]

    if adata_rna_DT.n_obs == 0 or adata_met_DT.n_obs == 0:
        raise ValueError("One of the datasets has no valid samples after filtering.")


    batch_size = min(50, adata_rna_DT.n_obs, adata_met_DT.n_obs)
    # OT MET->RNA
    adata_rna_DT_OT = adata_rna_DT
    print("Analyzing dataset")
    adata_rna_DT_OT = up.Run(adatas=[adata_rna_DT, adata_met_DT], mode='v', iteration=500, batch_size=batch_size)

    # # OT RNA->MET
    adata_mix = adata_rna_DT_OT[adata_rna_DT_OT.obs['type'] == 'sample']
    adata_ref = adata_rna_DT_OT[adata_rna_DT_OT.obs['type'] == 'ref']

    # Convert each subset to a pandas DataFrame
    def adata_to_df(adata):
        """Convert AnnData object to dense pandas DataFrame"""
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
        return pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)

    # Create the dictionary
    unit = {
        "mix": adata_to_df(adata_mix).T,
        "ref": adata_to_df(adata_ref).T,
        "ref_scRNA" : rna_unit["ref_scRNA"]
    }

    return unit
