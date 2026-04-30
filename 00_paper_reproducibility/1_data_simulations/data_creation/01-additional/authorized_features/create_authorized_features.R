# datasets = source("datasets.R")$value

datasets <- read.csv("datasets.csv", header = FALSE, comment.char = "#", col.names = c("datasets", "names_dataset"))[["datasets"]]


mixes = lapply(datasets, function(x) {
  print(x)
  readRDS(paste0("01_mixes/mixes1_",x,"_pdac.rds"))
})
names(mixes) = datasets

ref = readRDS("01_references/reference_pdac.rds")

authorized_genes = rownames(ref$ref_bulkRNA)
authorized_probes = rownames(ref$ref_met)
for (i in seq_along(datasets)) {
  authorized_genes = intersect(authorized_genes, rownames(mixes[[i]]$mix_rna))
  authorized_probes = intersect(authorized_probes, rownames(mixes[[i]]$mix_met))
}

length(authorized_genes)
length(authorized_probes)
saveRDS(authorized_genes, "01-additional/authorized_features/authorized_genes.rds")
saveRDS(authorized_probes, "01-additional/authorized_features/authorized_probes.rds")


library(IlluminaHumanMethylation27kanno.ilmn12.hg19)
pf27k = data.frame(getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19))
authorized_probes27k = intersect(authorized_probes, rownames(pf27k))
length(authorized_probes27k)
saveRDS(authorized_probes27k, "01-additional/authorized_features/authorized_probes27k.rds")


# for (i in seq_along(datasets)) {
#   tmp = mixes[[i]]
#   tmp$mix_rna = tmp$mix_rna[authorized_genes,]
#   tmp$mix_met = tmp$mix_met[authorized_probes,]
#   saveRDS(tmp, file = paste0("authorized_features/mixes1_",datasets[i],"_filtered_pdac.rds"))
# }

# tmp = ref
# tmp$ref_bulkRNA = tmp$ref_bulkRNA[authorized_genes,]
# tmp$ref_met = tmp$ref_met[authorized_probes,]
# saveRDS(tmp, file = "reference_filtered_pdac.rds")
