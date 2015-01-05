# code to combine p values from enrichment and proximity clustering of de novos

#' function to combine p values
#' 
#' @param x vector of P values for a gene
#' @export
#' 
#' @return single P value for gene
fishersMethod <- function(x) {
    x = x[!is.na(x)]
    if (length(x) == 0) { return(NA) }
    
    adjusted = pchisq(-2 * sum(log(x)), df=2 * length(x), lower=FALSE)
    
    return(adjusted)
}

#' combine P values from enrichment and clustering tests into a single P value
#' 
#' @param enrichment_path path to enrichment testing results
#' @param enrichment_path path to clustering testing results
#' @export
#' 
#' @return a merged dataset where the P values have been combined
combine_analyses <- function(enrichment_path, clustering_path) {
    # read in p values from clustering analysis, only for genes with >1 mutation
    clust = read.table(clustering_path, header=TRUE, sep="\t")
    clust = cast(clust, gene_id ~ mutation_category, value="probability", mean)
    
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
    
    # adjust the P values by false-discovery rate
    num.tests = 18500
    merged$combined_fdr = p.adjust(merged$p_combined, method="BH", n=num.tests)
    
    return(merged)
}

#' example of aanalysis using the above functions 
main <- function() {
    CODE_DIR = "/nfs/users/nfs_j/jm33/apps/mupit"
    RESULTS_DIR = file.path(CODE_DIR, "results")
    clustering_path = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.clustering_results.txt")
    enrichment_path = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.enrichment_results.txt")
    output_path = file.path(RESULTS_DIR, "de_novos.ddd_4k.meta-analysis.combined.txt")
    
    combine_analyses(enrichment_path, clustering_path, output_path)
    write.table(merged, file=output_path, row.names=FALSE, quote=FALSE, sep="\t")
}


