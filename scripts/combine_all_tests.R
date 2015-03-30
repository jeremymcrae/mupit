# code to combine the different tests statistics together into a single table

library(mupit)

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/mupit"
RESULTS_DIR = file.path(CODE_DIR, "results")

DDG2P_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/DDG2P_freeze_with_gencode19_genomic_coordinates_20141118_fixed.txt"

include_ddg2_status <- function(merged) {
    # load the DDG2P, so that we can annotated the genes with their DDG2P status
    ddg2p = read.table(DDG2P_PATH, sep="\t", header=TRUE, stringsAsFactor=FALSE)
    ddg2p = ddg2p[ddg2p$DDD_category != "Possible DD Gene", ]
    
    # annotate each column with DDG2P status
    merged$in_ddg2p = merged$hgnc %in% ddg2p$gencode_gene_name
    
    return(merged)
}

diagnosed = c("with_diagnosed", "without_diagnosed")
for (diagnosed_string in diagnosed) {
    prefix = paste("de_novos.ddd_4k.", diagnosed_string, sep="")

    # define the normal de novo analysis result files
    meta_clust = file.path(RESULTS_DIR, paste(prefix, ".meta-analysis.clustering_results.txt", sep=""))
    meta_enrich = file.path(RESULTS_DIR, paste(prefix, ".meta-analysis.enrichment_results.txt", sep=""))
    clust = file.path(RESULTS_DIR, paste(prefix, ".ddd_only.clustering_results.txt", sep=""))
    enrich = file.path(RESULTS_DIR, paste(prefix, ".ddd_only.enrichment_results.txt", sep=""))
    pheno = file.path(RESULTS_DIR, paste(prefix, ".ddd_only.phenotype_similarity_results.txt", sep=""))

    merged = combine_tests(meta_clust, meta_enrich, clust, enrich, pheno)
    merged = include_ddg2_status(merged)
    write.table(merged, file=file.path(RESULTS_DIR, paste(prefix, ".all.txt", sep="")), row.names=F, quote=F, sep="\t")
}
