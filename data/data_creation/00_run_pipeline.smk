
def c(*args): return list(args)

NULL = None



import csv
datasets = {row[0].strip(): row[1].strip() for row in csv.reader(open("datasets.csv")) if not row[0].startswith("#")}.keys()


reference_files = ["01_references/reference_pdac.rds"]
mixes_files = [f"01_mixes/mixes{dc}_{dataset}_pdac.rds" for dataset in datasets for dc in [1,2]]
filteredmixes_files = [f"01_mixes/filteredmixes{dc}_{dataset}_pdac.rds" for dataset in datasets for dc in [1,2]]
filteredreference_files = ["01_references/filteredreference_pdac.rds"]
groundtruth_files = [f"01_groundtruth/groundtruth{dc}_{dataset}_pdac.rds" for dataset in datasets for dc in [1,2]]


rule target:
    threads: 1
    message: "-- Rule target completed --"
    input:
      reference_files,
      mixes_files,
      groundtruth_files,
      filteredmixes_files,
      filteredreference_files

rule datasets:
    threads: 1
    message: "-- Rule target completed --"
    input:
      reference_files,
      mixes_files,
      groundtruth_files,
      filteredmixes_files,
      filteredreference_files


rule source:
    input: 
      rmd_file = "01.1_src_pdac.Rmd",
    output: 
      source_file = "01_references/source_pdac.rds",
    threads: 1
    log: file = "logs/01.1_src_pdac.Rout"
    shell:"""
mkdir -p tmp_all
mkdir -p 01_html
mkdir -p 01_references
mkdir -p 01_groundtruth
mkdir -p 01_mixes
RCODE="source_file='{output.source_file}'; rmarkdown::render('{input.rmd_file}', output_file='01.1_src_pdac.html');"
echo $RCODE | Rscript - 2>&1 > {log.file}
mv 01.1_src_pdac.html 01_html/01.1_src_pdac.html
"""

rule reference:
    input: 
      rmd_file = "01.1_ref_pdac.Rmd",
    output: 
      reference_file = "01_references/reference_pdac.rds",
    threads: 1
    log: file = "logs/01.1_ref_pdac.Rout"
    shell:"""
mkdir -p tmp_all
mkdir -p 01_html
mkdir -p 01_references
mkdir -p 01_groundtruth
mkdir -p 01_mixes
RCODE="reference_file='{output.reference_file}'; rmarkdown::render('{input.rmd_file}', output_file='01.1_ref_pdac.html');"
echo $RCODE | Rscript - 2>&1 > {log.file}
mv 01.1_ref_pdac.html 01_html/01.1_ref_pdac.html
"""

rule mixes_and_groundtruth:
    input: 
      rmd_file = "01.2_dataset/01.2_dataset_{dataset}.Rmd",
      source_file = "01_references/source_{domain}.rds",
    output: 
      mixes_file1 = "01_mixes/mixes1_{dataset}_{domain}.rds",
      mixes_file2 = "01_mixes/mixes2_{dataset}_{domain}.rds",
      groundtruth_file1 = "01_groundtruth/groundtruth1_{dataset}_{domain}.rds",
      groundtruth_file2 = "01_groundtruth/groundtruth2_{dataset}_{domain}.rds",
    threads: 1
    log: file = "logs/01.2_dataset_{dataset}_{domain}.Rout"
    shell:"""
RCODE="source_file='{input.source_file}'; mixes_files=c('{output.mixes_file1}','{output.mixes_file2}'); groundtruth_files=c('{output.groundtruth_file1}','{output.groundtruth_file2}'); rmarkdown::render('{input.rmd_file}', output_file='01.2_dataset_{wildcards.dataset}_{wildcards.domain}.html', intermediates_dir='tmp_all/{wildcards.dataset}_{wildcards.domain}',clean = FALSE);"
echo $RCODE | Rscript - 2>&1 > {log.file}
mv 01.2_dataset/01.2_dataset_{wildcards.dataset}_{wildcards.domain}.html 01_html/01.2_dataset_{wildcards.dataset}_{wildcards.domain}.html
"""
# Warning source_file is called reference_file in all 01.2_dataset/01.2_dataset_XXX.Rmd

rule filteredmixes:
    input: 
      rmd_file = "01.2_filterdataset.Rmd",
      mixes_file1 = "01_mixes/mixes1_{dataset}_{domain}.rds",
      mixes_file2 = "01_mixes/mixes2_{dataset}_{domain}.rds",
      authorized_probes_rds = "01-additional/authorized_features/authorized_probes27k.rds",
      authorized_genes_rds = "01-additional/authorized_features/authorized_genes.rds",
    output: 
      mixes_file1 = "01_mixes/filteredmixes1_{dataset}_{domain}.rds",
      mixes_file2 = "01_mixes/filteredmixes2_{dataset}_{domain}.rds",
    threads: 1
    log: file = "logs/01.2_filterdataset_{dataset}_{domain}.Rout"
    shell:"""
RCODE="domain='{wildcards.domain}'; authorized_genes_file='{input.authorized_genes_rds}' ; authorized_probes_file='{input.authorized_probes_rds}' ; mixes_files=c('{input.mixes_file1}','{input.mixes_file2}'); rmarkdown::render('{input.rmd_file}', output_file='01.2_filterdataset_{wildcards.dataset}_{wildcards.domain}.html', intermediates_dir='tmp_all/{wildcards.domain}',clean = FALSE);"
echo $RCODE | Rscript - 2>&1 > {log.file}
mv 01.2_filterdataset_{wildcards.dataset}_{wildcards.domain}.html 01_html/01.2_filterdataset_{wildcards.dataset}_{wildcards.domain}.html
"""

rule filteredreference:
    input:
      rmd_file = "01.2_filterreference.Rmd",
      reference_file = "01_references/reference_{domain}.rds",
      authorized_probes_rds = "01-additional/authorized_features/authorized_probes27k.rds",
      authorized_genes_rds = "01-additional/authorized_features/authorized_genes.rds"
    output:
      filteredreference_file = "01_references/filteredreference_{domain}.rds"
    threads: 1
    log: file = "logs/01.2_filterreference_{domain}.Rout"
    shell:"""
RCODE="domain='{wildcards.domain}'; authorized_genes_file='{input.authorized_genes_rds}' ; authorized_probes_file='{input.authorized_probes_rds}' ;reference_file='{input.reference_file}'; rmarkdown::render('{input.rmd_file}', output_file='01.2_filterreference_{wildcards.domain}.html', intermediates_dir='tmp_all/{wildcards.domain}',clean = FALSE);"
echo $RCODE | Rscript - 2>&1 > {log.file}
mv 01.2_filterreference_{wildcards.domain}.html 01_html/01.2_filterreference_{wildcards.domain}.html
"""

rule clean:
    threads: 1
    shell:"""
rm -rf logs/
rm -rf tmp_all/
rm -rf 01_*
rm -rf 02_prediction/
rm -rf 03_scores/
rm -rf 04_visu/
rm -rf datashare/
"""

rule gantt:
    threads: 1
    shell:"""
smgantt
"""