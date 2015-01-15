# code to combine p values from enrichment and proximity clustering of de novos,
# as well as to get the most significant P value from different subsets of tests

#' function to combine p values, using Fisher's method
#' 
#' @param x vector of P values for a gene
#' @export
#' 
#' @return single P value for gene
fishersMethod <- function(x) {
    x = x[!is.na(x)]
    if (length(x) == 0) { return(NA) }
    
    adjusted = pchisq(-2 * sum(log(x)), df=2 * length(x), lower.tail=FALSE)
    
    return(adjusted)
}

#' combine P values from enrichment and clustering tests into a single P value
#' 
#' @param enrichment_path path to enrichment testing results
#' @param clustering_path path to clustering testing results
#' @export
#' 
#' @return a merged dataset where the P values have been combined
combine_enrichment_and_clustering <- function(enrichment_path, clustering_path) {
    # read in p values from clustering analysis, only for genes with >1 mutation
    clust = read.table(clustering_path, header=TRUE, sep="\t")
    clust = reshape::cast(clust, gene_id ~ mutation_category, value="probability", mean)
    
    # read in p values from mupit analyses
    enriched = read.table(enrichment_path, header=TRUE, sep="\t")
    
    # merge the datasets
    merged = merge(enriched, clust, by.x="hgnc", by.y="gene_id", all.x=TRUE)
    
    # calculate a combined p-value for each gene, using the synonymous and
    # missense categories from each of the enrichment and clustering analyses.
    # We don't expect the loss-of-function de novos to be clustered, so we don't
    # use those. This way we can later get the minimum of the loss-of-function
    # tests, and the minimum of the combined P values.
    p_values = merged[, c("synonymous", "missense", "p_func", "p_synonymous")]
    merged$p_combined = apply(p_values, 1, fishersMethod)
    
    # # we might want to state the correlation between the P values from the
    # # different analyses, to show that the two analyses are independent.
    # print(cor(-log10(p_values$synonymous), -log10(p_values$p_synonymous), use="complete.obs")^2)
    # print(cor(-log10(p_values$missense), -log10(p_values$p_func), use="complete.obs")^2)
    
    # adjust the P values by false-discovery rate
    num.tests = 18500
    merged$combined_fdr = p.adjust(merged$p_combined, method="BH", n=num.tests)
    
    return(merged)
}

#' find the most significant P value for each gene from the P values from 
#' different subsets and different tests
#' 
#' @param meta_clust path to clustering results for the meta-analysis subset
#' @param meta_enrich path to enrichment results for the meta-analysis subset
#' @param clust path to clustering results for the ddd only subset
#' @param enrich path to enrichment results for the ddd only subset
#' @export
#' 
#' @return data frame with the columns from all the datasets, as well as minimum
#'     P values from each subset for each gene, and overall minimum P values for
#'     each gene.
combine_tests <- function(meta_clust, meta_enrich, clust, enrich) {
    # load all the data files in. Previously Matt was checking two different 
    # subsets, one where the undiagnosed individuals were included, and one  
    # where they were excluded. I've swapped to only using the undiagnosed 
    # excluded subset
    meta = combine_enrichment_and_clustering(meta_enrich, meta_clust) 
    ddd = combine_enrichment_and_clustering(enrich, clust) 
    
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

#' example of analysis using the above functions 
main <- function() {
    CODE_DIR = "/nfs/users/nfs_j/jm33/apps/mupit"
    RESULTS_DIR = file.path(CODE_DIR, "results")
    clustering_path = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.clustering_results.txt")
    enrichment_path = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.enrichment_results.txt")
    output_path = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.combined.txt")
    
    merged = combine_enrichment_and_clustering(enrichment_path, clustering_path)
    write.table(merged, file=output_path, row.names=FALSE, quote=FALSE, sep="\t")
}


