# code to combine the different tests statistics together into a single table

library(mupit)

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/mupit"
RESULTS_DIR = file.path(CODE_DIR, "results")

diagnosed_string = "without_diagnosed"
prefix = paste("de_novos.ddd_4k.", diagnosed_string, sep="")

# define the normal de novo analysis result files
ddd_meta_clust = file.path(RESULTS_DIR, paste(prefix, ".meta-analysis.clustering_results.txt", sep=""))
ddd_meta_enrich = file.path(RESULTS_DIR, paste(prefix, ".meta-analysis.enrichment_results.txt", sep=""))
ddd_clust = file.path(RESULTS_DIR, paste(prefix, ".ddd_only.clustering_results.txt", sep=""))
ddd_enrich = file.path(RESULTS_DIR, paste(prefix, ".ddd_only.enrichment_results.txt", sep=""))

merged = combine_tests(ddd_meta_clust, ddd_meta_enrich, ddd_clust, ddd_enrich)
write.table(merged, file=file.path(RESULTS_DIR, paste(prefix, ".all.txt", sep="")), row.names=F, quote=F, sep="\t")
