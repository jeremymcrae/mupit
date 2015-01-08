# code to combine the different tests statistics together into a single table

library(mupit)

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/mupit"
RESULTS_DIR = file.path(CODE_DIR, "results")

#' find the most significant P value for each gene from the P values from 
#' different subsets and different tests
#' 
#' @param meta_clust path to clustering results for the meta-analysis subset
#' @param meta_enrich path to enrichment results for the meta-analysis subset
#' @param clust path to clustering results for the ddd only subset
#' @param enrich path to enrichment results for the ddd only subset
#' 
#' @return data frame with the columns from all the datasets, as well as minimum
#'     P values from each subset for each gene, and overall minimum P values for
#'     each gene.
get_most_significant <- function(meta_clust, meta_enrich, clust, enrich) {
    # load all the data files in. Previously Matt was checking two different 
    # subsets, one where the undiagnosed individuals were included, and one  
    # where they were excluded. I've swapped to only using the undiagnosed 
    # excluded subset
    meta = combine_analyses(meta_enrich, meta_clust) 
    ddd = combine_analyses(enrich, clust) 
    
    # need names that are more informative as same across files, add prefix
    names(ddd) = paste("ddd", names(ddd), sep=".")
    names(meta) = paste("meta", names(meta), sep=".")
    
    # merge together files, focusing on genes with DNMs in DDD
    merged = merge(meta, ddd, by.x="meta.hgnc", by.y="ddd.hgnc", all.x=TRUE)
    
    # calculate minimum p value across LoF and func + clustering tests for each dataset
    merged$p_min_ddd = apply(merged[, c("ddd.p_lof", "ddd.p_combined")], 1, min, na.rm=TRUE)
    merged$p_min_meta = apply(merged[, c("meta.p_lof",  "meta.p_combined")], 1, min, na.rm=TRUE)
    merged$p_min = apply(merged[, c("p_min_ddd", "p_min_meta")], 1, min, na.rm=TRUE)
    
    return(merged)
}

# define the normal de novo analysis result files
ddd_meta_clust = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.clustering_results.txt")
ddd_meta_enrich = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.enrichment_results.txt")
ddd_clust = file.path(RESULTS_DIR, "de_novos.ddd_4k.ddd_only.enrichment_results.txt")
ddd_enrich = file.path(RESULTS_DIR, "de_novos.ddd_4k.ddd_only.clustering_results.txt")

# define the seizure analysis result files
seizures_meta_clust = file.path(RESULTS_DIR, "de_novos.seizures.meta-analysis.clustering_results.txt")
seizures_meta_enrich = file.path(RESULTS_DIR, "de_novos.seizures.meta-analysis.enrichment_results.txt")
seizures_clust = file.path(RESULTS_DIR, "de_novos.seizures.ddd_only.clustering_results.txt")
seizures_enrich = file.path(RESULTS_DIR, "de_novos.seizures.ddd_only.enrichment_results.txt")


meta_clust = ddd_meta_clust
meta_enrich = ddd_meta_enrich
clust = ddd_clust
enrich = ddd_enrich

merged = get_most_significant(ddd_meta_clust, ddd_meta_enrich, ddd_clust, ddd_enrich)
write.table(merged, file=file.path(RESULTS_DIR, "de_novos.ddd_4k.all.txt"), row.names=F, quote=F, sep="\t")


merged = get_most_significant(seizures_meta_clust, seizures_meta_enrich, seizures_clust, seizures_enrich)
write.table(merged, file=file.path(RESULTS_DIR, "de_novos.seizures.all.txt"), row.names=F, quote=F, sep="\t")


