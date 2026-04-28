# Create the list of association between methylation sites and genes

probes_feature = readRDS("01-additional/probes_features.rds")

# withdraw all CpG sites in intergenic region
probes_feature = probes_feature[probes_feature$gene != "", ]
