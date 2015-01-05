# code to combine the different tests statistics together into a single table

library(mupit)

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/mupit"
RESULTS_DIR = file.path(CODE_DIR, "results")

# load all the data files in. Previously Matt was checking two different 
# subsets, one where the undiagnosed individuals were included, and one where 
# they were excluded. I've swapped to only using the undiagnosed excluded subset
meta_clust = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.clustering_results.txt")
meta_enrich = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.enrichment_results.txt")
meta.undiagnosed = combine_analyses(meta_enrich, meta_clust) 

ddd_clust = file.path(RESULTS_DIR, "de_novos.ddd_4k.ddd_only.clustering_results.txt")
ddd_enrich = file.path(RESULTS_DIR, "de_novos.ddd_4k.ddd_only.enrichment_results.txt")
ddd.undiagnosed = combine_analyses(meta_enrich, meta_clust) 

# need names that are more informative as same across files, add prefix
# names(ddd.all) = paste("ddd.all", names(ddd.all), sep=".")
# names(meta.all) = paste("meta.all", names(meta.all), sep=".")
names(ddd.undiagnosed) = paste("ddd.undiagnosed", names(ddd.undiagnosed), sep=".")
names(meta.undiagnosed) = paste("meta.undiagnosed", names(meta.undiagnosed), sep=".")

# merge together files, focusing on genes with DNMs in DDD
merged = merge(meta.undiagnosed, ddd.undiagnosed, by.x="hgnc", by.y="hgnc", all.x=TRUE)

# calculate minimum p value across LoF and func + clustering tests for each dataset
min_ddd_undiagnosed_p = apply(merged[, c("ddd.undiagnosed.p_lof", "ddd.undiagnosed.p_combined")], 1, min, na.rm=TRUE)
min_meta_undiagnosed_p = apply(merged[, c("meta.undiagnosed.p_lof",  "meta.undiagnosed.p_combined")], 1, min, na.rm=TRUE)
min_undiagnosed_p = apply(cbind(min_ddd_undiagnosed_p, min_meta_undiagnosed_p), 1, min, na.rm=TRUE)

# write out full table
merged = cbind(merged, min_ddd_undiagnosed_p, min_meta_undiagnosed_p, min_undiagnosed_p)
write.table(merged, file="all.tests.combined.txt", row.names=F, quote=F, sep="\t")
