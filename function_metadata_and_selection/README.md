# Functions metadata


This folder contains various setups for selecting and describing functions to be included in the Nextflow pipeline.
Folder contents:

    * `CI_functions_selection`: Used for running the pipeline with the resume CI.
    * `limited_EI`: Created to test with a subset of the full dataset, generating all combinations. All selections (except dataset.yml) are symbolic links from the `template_complete_run` folder.
    * `setup`: Contains several setups used for benchmarking Snakemake vs. Python, or for benchmarking specific functions.
    * `Local_test_setup`: A small subset of function selections used for the Full CI and local pipeline testing.


## Benchamark to test snakemake vs nextflow 
/!\ Change branch to `smk_nf_used_for_benchmark`

```
cd function_metadata_and_selection
conda activate hadaca3framework_env
python benchmark.py
```

## Run the analysis: 

```
cd function_metadata_and_selection
conda activate hadaca3framework_env
Rscript -e "rmarkdown::render('analysis.Rmd')"
```

## setup details : 

### setup1 
```
job                             count
----------------------------  -------
cleaning_mix                        1
cleaning_ref                        1
features_selection                  5
late_integration                    4
metaanalysis                        1
prediction_deconvolution_met        2
prediction_deconvolution_rna        2
preprocessing                       5
scoring                             4
total                              25
```


### setup2 
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                  7
late_integration                    8
metaanalysis                        1
prediction_deconvolution_met        4
prediction_deconvolution_rna        4
preprocessing                       7
scoring                             8
total                              42
```

### setup3
```
job                             count
----------------------------  -------
cleaning_mix                        4
cleaning_ref                        1
features_selection                 11
late_integration                   16
metaanalysis                        1
prediction_deconvolution_met        8
prediction_deconvolution_rna        8
preprocessing                      11
scoring                            16
total                              76

```

### setup4 
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 16
late_integration                  128
metaanalysis                        1
prediction_deconvolution_met       64
prediction_deconvolution_rna        4
preprocessing                      10
scoring                           128
total                             354

```

### setup5 

```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 16
late_integration                  384
metaanalysis                        1
prediction_deconvolution_met       64
prediction_deconvolution_rna        4
preprocessing                      10
scoring                           384
total                             866

```

### setup6
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 16
late_integration                  512
metaanalysis                        1
prediction_deconvolution_met       64
prediction_deconvolution_rna        4
preprocessing                      10
scoring                           512
total                            1122
```



### setup7
(4min45 dag creation dry run )
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 26
late_integration                 4608
metaanalysis                        1
prediction_deconvolution_met       96
prediction_deconvolution_rna       96
preprocessing                      13
scoring                          4608
total                            9451
```


### setup8
(13min dag creation dry run)
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 33
late_integration                10368
metaanalysis                        1
prediction_deconvolution_met      216
prediction_deconvolution_rna       96
preprocessing                      13
scoring                         10368
total                           21098
```

### setup9
(43min dag creation)
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 33
late_integration                31104
metaanalysis                        1
prediction_deconvolution_met      216
prediction_deconvolution_rna       96
preprocessing                      13
scoring                         31104
total                           62570
```


### setup 10 

2,5, 4,2 6 for the number datasets, pp, fs, de, li functions respectively. 

(85min dag creation)
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 33
late_integration                62208
metaanalysis                        1
prediction_deconvolution_met      216
prediction_deconvolution_rna       96
preprocessing                      13
scoring                         62208
total                          124778
```

## setup 11 (not benchmarked)
Only one li fonction ! 

```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 33
late_integration                10368
metaanalysis                        1
prediction_deconvolution_met      216
prediction_deconvolution_rna       96
preprocessing                      13
scoring                         10368
total                           21098
```

## setup all_any 
nb de fonction : 
datasets 2 
pp : 3
fs  :2
dec : 2
li : 1

(125m37,285s)

```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 42
late_integration               62 208
metaanalysis                        1
prediction_deconvolution_met      144
prediction_deconvolution_rna      864
preprocessing                      21
scoring                        62 208
total                         125 491
```
