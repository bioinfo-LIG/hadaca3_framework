# hadaca3 data creation 


## Conda environement

Set up your conda environement as follow:

```
conda create -n hadaca3data_env
conda activate hadaca3data_env

mamba install -y  -c bioconda -c conda-forge -c r r-base snakemake=7.32.4 python=3.9 r-rmarkdown r-seurat r-clue r-ggpubr  r-viridis   

```

graphviz python-kaleido tenacity plotly r-bisquerna r-extraDistr r-MASS r-EPIC r-fmsb bioconductor-toast bioconductor-omicade4 r-mixomics r-mixkernel rpy2 scikit-learn keras tensorflow bioconductor-viper bioconductor-ADImpute r-WGCNA r-see r-ggfortify scanpy bioconda::bioconductor-graph  bioconda::r-wgcna  bioconda::bioconductor-edger bioconda::bioconductor-rgraphviz



## How to start?

To execute the entire pipeline (generate sources, references, datsets and run baselines, scoring and visualisation) execute the following lines. You will probably need to get original data, so refer to the next section if needed.

```
snakemake --cores 1 -s 00_run_pipeline.py -p clean  # remove .rds files
snakemake --cores 4 -s 00_run_pipeline.py -pn       # dry-run
```

This pipeline can be visualised by generating its DAG:

```
snakemake --forceall --dag -s 00_run_pipeline.py | dot -Tpdf > dag.pdf
```


## Configure your work. env. by getting original data

The section describes which data are needed to execute the entire pipeline and provide the code to download it.


List of files needed for the correct execution : 
* ~/projects/pancreas_deconv_data/results/20210624/02.1_output/peng_pool_labelled.rds
* ~/projects/datashare/cometh_lot1/*  (mag)
* ~/projects/datashare/data_cometh_lot2/*  (mag)
* ~/projects/datashare/input_ref_hadaca3/* (lucie)
* ~/projects/datashare/genref_hadaca3/genref_hadaca2_rna.rds
* ~/projects/datashare/genref_hadaca3/genref_hadaca2_meth.rds
* ~/projects/datashare/genref_hadaca3/00_peng_k_2019.rds
* ~/projects/datashare/genref_hadaca3/00_raghavan_s_2021.rds
* ~/projects/datashare/genref_hadaca3/00_baron_m_2016.rds
* ~/projects/datashare/genref_hadaca3/cpg_probes_400k
```
scp -r dahu.ciment:/bettik/lamothlu/projects/datashare/input_ref_hadaca3 ~/projects/datasahre/
mkdir -p ~/projects/datashare/cometh_lot1/transcriptome
rsync -auvP dahu.ciment:/home/richamag/projects/datashare/cometh_lot1/transcriptome/gene_counts_normalised.txt ~/projects/datashare/cometh_lot1/transcriptome/gene_counts_normalised.txt

mkdir -p ~/projects/datashare/cometh_lot1/methylation
rsync -auvP dahu.ciment:/home/richamag/projects/datashare/cometh_lot1/methylation/test_data_met_new.rds ~/projects/datashare/cometh_lot1/methylation/test_data_met_new.rds

rsync -auvP dahu.ciment:/home/richamag/projects/datashare/cometh_lot1/test_solution.rds ~/projects/datashare/cometh_lot1/test_solution.rds

rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/input_ref_hadaca3/ ~/projects/datashare/input_ref_hadaca3


mkdir -p ~/projects/datashare/cometh_lot2/methylation/
mkdir -p ~/projects/datashare/cometh_lot2/transcriptome/
mkdir -p ~/projects/datashare/cometh_lot2/prediction_owkin/
rsync -auvP dahu.ciment:/bettik/amblaeli/projects/datashare/cometh_lot2/methylation/test_data_met_new.rds          ~/projects/datashare/cometh_lot2/methylation/     
rsync -auvP dahu.ciment:/bettik/amblaeli/projects/datashare/cometh_lot2/transcriptome/gene_counts_normalised.txt   ~/projects/datashare/cometh_lot2/transcriptome/   
rsync -auvP dahu.ciment:/bettik/amblaeli/projects/datashare/cometh_lot2/prediction_owkin/preds_clean_50_add_mr.csv ~/projects/datashare/cometh_lot2/prediction_owkin/

rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/pseudobulk_mix_hadaca3/ ~/projects/datashare/pseudobulk_mix_hadaca3

rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/genref_hadaca3/00_peng_k_2019.rds        ~/projects/datashare/genref_hadaca3/
rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/genref_hadaca3/00_raghavan_s_2021.rds    ~/projects/datashare/genref_hadaca3/
rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/genref_hadaca3/00_baron_m_2016.rds       ~/projects/datashare/genref_hadaca3/
rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/genref_hadaca3/genref_hadaca2_meth.rds   ~/projects/datashare/genref_hadaca3/
rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/genref_hadaca3/genref_hadaca2_rna.rds    ~/projects/datashare/genref_hadaca3/
rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/genref_hadaca3/00_peng_k_2019.rds        ~/projects/datashare/genref_hadaca3/
rsync -auvP dahu.ciment:/bettik/richamag/projects/pancreas_deconv_data/results/20210624/02.1_output/peng_pool_labelled.rds ~/projects/datashare/pancreas_deconv_data/results/20210624/02.1_output/

rsync -auvP dahu.ciment:/bettik/lamothlu/projects/datashare/genref_hadaca3/lin_pool_labelled.rds  ~/projects/datashare/genref_hadaca3/
```

01.2_dataset_insilicopseudobulk.Rmd takes took long to run (~ 2h) 
For now you will need to rsync : 

~/projects/datashare/genref_hadaca3/pseudobulk_mix/* (Lucie) 

If you want to run all of the chunks of the code you will still need to rsync this file for pseudo sc meth generation:

~/projects/datashare/genref_hadaca3/cpg_probes_400k


Theses files or folder can be retrieved for instance with the command : 
```
```