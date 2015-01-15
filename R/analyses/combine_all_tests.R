# code to combine the different tests statistics together into a single table

library(mupit)

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/mupit"
RESULTS_DIR = file.path(CODE_DIR, "results")

# define the normal de novo analysis result files
ddd_meta_clust = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.clustering_results.txt")
ddd_meta_enrich = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.enrichment_results.txt")
ddd_clust = file.path(RESULTS_DIR, "de_novos.ddd_4k.ddd_only.clustering_results.txt")
ddd_enrich = file.path(RESULTS_DIR, "de_novos.ddd_4k.ddd_only.enrichment_results.txt")

merged = combine_tests(ddd_meta_clust, ddd_meta_enrich, ddd_clust, ddd_enrich)
write.table(merged, file=file.path(RESULTS_DIR, "de_novos.ddd_4k.all.txt"), row.names=F, quote=F, sep="\t")
