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
    names(clust)[2:3] = paste("p_", names(clust)[2:3], "_clust", sep="")
    
    # read in p values from mupit analyses
    enriched = read.table(enrichment_path, header=TRUE, sep="\t")
    
    # merge the datasets
    merged = merge(enriched, clust, by.x="hgnc", by.y="gene_id", all.x=TRUE)
    
    # calculate a combined p-value for each gene, using the synonymous and
    # missense categories from each of the enrichment and clustering analyses.
    # We don't expect the loss-of-function de novos to be clustered, so we don't
    # use those. This way we can later get the minimum of the loss-of-function
    # tests, and the minimum of the combined P values.
    p_values = merged[, c("p_missense_clust", "p_func")]
    merged$p_combined = apply(p_values, 1, fishersMethod)
    
    # # we might want to state the correlation between the P values from the
    # # different analyses, to show that the two analyses are independent.
    # print(cor(-log10(p_values$missense), -log10(p_values$p_func), use="complete.obs")^2)
    
    # adjust the P values by false-discovery rate
    num.tests = 18500
    merged$combined_fdr = p.adjust(merged$p_combined, method="BH", n=num.tests)
    
    # calculate minimum p value across LoF and func + clustering tests
    merged$p_min = apply(merged[, c("p_lof", "p_combined")], 1,
        function(x) ifelse(length(x[!is.na(x)]) > 0, min(x[!is.na(x)]), NA))
    
    return(merged)
}

#' find the most significant P value for each gene from the P values from
#' different subsets and different tests
#'
#' @param meta_clust path to clustering results for the meta-analysis subset
#' @param meta_enrich path to enrichment results for the meta-analysis subset
#' @param clust path to clustering results for the ddd only subset
#' @param enrich path to enrichment results for the ddd only subset
#' @param phenotype path to phenotype similarity testing results
#' @export
#'
#' @return data frame with the columns from all the datasets, as well as minimum
#'     P values from each subset for each gene, and overall minimum P values for
#'     each gene.
combine_tests <- function(meta_clust, meta_enrich, clust, enrich, phenotype) {
    # load all the data files in
    meta = combine_enrichment_and_clustering(meta_enrich, meta_clust)
    ddd = combine_enrichment_and_clustering(enrich, clust)
    
    # read in p values from HPO similarity analyses
    phenotypes = read.table(phenotype, header=TRUE, sep="\t")
    
    ddd = merge(ddd, phenotypes, by="hgnc", all.x=TRUE)
    ddd$p_min_with_pheno = apply(ddd[, c("p_min", "hpo_similarity_p_value")], 1, fishersMethod)
    
    # need names that are more informative as same across files, add prefix
    names(ddd)[3:length(names(ddd))] = paste("ddd", names(ddd)[3:length(names(ddd))], sep=".")
    names(meta)[3:length(names(meta))] = paste("meta", names(meta)[3:length(names(meta))], sep=".")
    
    # merge together files, focusing on genes with DNMs in DDD
    merged = merge(meta, ddd, by=c("hgnc", "chrom"), all.x=TRUE)
    merged$p_min = apply(merged[, c("ddd.p_min_with_pheno", "meta.p_min")], 1,
        function(x) ifelse(length(x[!is.na(x)]) > 0, min(x[!is.na(x)]), NA))
    
    return(merged)
}
