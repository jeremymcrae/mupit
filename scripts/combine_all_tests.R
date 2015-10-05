# code to combine the different tests statistics together into a single table

library(mupit)

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/mupit"
RESULTS_DIR = file.path(CODE_DIR, "results")

DDG2P_PATH = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"

include_ddg2p_status <- function(merged, ddg2p_path) {
    # load the DDG2P, so that we can annotated the genes with their DDG2P status
    ddg2p = load_ddg2p(ddg2p_path)
    
    # annotate each column with DDG2P status
    merged$in_ddg2p = merged$hgnc %in% ddg2p$gene
    merged$in_ddg2p_dominant = merged$hgnc %in% ddg2p$gene[ddg2p$dominant]
    
    return(merged)
}

diagnosed = c("with_diagnosed", "without_diagnosed")
for (diagnosed_string in diagnosed) {
    prefix = paste("de_novos.ddd_4k.", diagnosed_string, sep="")

    # define the normal de novo analysis result files
    meta_clust = file.path(RESULTS_DIR, paste(prefix, ".meta-analysis.clustering_results.txt", sep=""))
    meta_enrich = file.path(RESULTS_DIR, paste(prefix, ".meta-analysis.enrichment_results.2015-10-02.txt", sep=""))
    clust = file.path(RESULTS_DIR, paste(prefix, ".ddd_only.clustering_results.txt", sep=""))
    enrich = file.path(RESULTS_DIR, paste(prefix, ".ddd_only.enrichment_results.2015-10-02.txt", sep=""))
    pheno = file.path(RESULTS_DIR, paste(prefix, ".ddd_only.phenotype_similarity_results.txt", sep=""))

    merged = combine_tests(meta_clust, meta_enrich, clust, enrich, pheno)
    merged = include_ddg2p_status(merged, DDG2P_PATH)
    write.table(merged, file=file.path(RESULTS_DIR, paste(prefix, "all", Sys.Date(), "txt", sep=".")), row.names=F, quote=F, sep="\t")
}
