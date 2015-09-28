# code to combine p values from enrichment and proximity clustering of de novos,
# as well as to get the most significant P value from different subsets of tests

#' function to combine p values, using Fisher's method
#'
#' @param x vector of P values for a gene
#' @export
#'
#' @return single P value for gene
#' @examples
#' fishersMethod(c(0.01, 0.001))
#' fishersMethod(c(0.01, NA))
#' fishersMethod(c(NA, NA))
#' fishersMethod(c(0.01))
fishersMethod <- function(x) {
    x = x[!is.na(x)]
    if (length(x) == 0) { return(NA) }
    
    adjusted = pchisq(-2 * sum(log(x)), df=2 * length(x), lower.tail=FALSE)
    
    return(adjusted)
}

#' combine P values from enrichment and clustering tests into a single P value
#'
#' @param enriched dataframe of de novo enrichment results
#' @param clust dataframe of de novo clustering results
#' @param num_tests number of tests to correct for multiple testing
#' @export
#'
#' @return a merged dataset where the P values have been combined
#' @examples
#' enriched = read.table(header=TRUE, text="
#'     hgnc  enrichment_p_value
#'     GENE1 0.01
#'     GENE2 0.0001")
#' clust = read.table(header=TRUE, text="
#'     gene_id mutation_category probability
#'     GENE1   missense          0.1
#'     GENE1   nonsense          0.1
#'     GENE2   missense          0.001
#'     GENE2   nonsense          0.001")
#' combine_enrichment_and_clustering(enriched, clust)
combine_enrichment_and_clustering <- function(enriched, clust, num_tests=18500) {
    # read in p values from clustering analysis, only for genes with >1 mutation
    clust = reshape::cast(clust, gene_id ~ mutation_category, value="probability", mean)
    names(clust)[2:3] = paste("p_", names(clust)[2:3], "_clust", sep="")
    
    # merge the datasets
    merged = merge(enriched, clust, by.x="hgnc", by.y="gene_id", all.x=TRUE)
    
    # calculate a combined p-value for each gene. We don't expect the
    # loss-of-function de novos to be clustered, so we don't use that.
    p_values = merged[, c("enrichment_p_value", "p_missense_clust")]
    merged$p_combined = apply(p_values, 1, fishersMethod)
    
    # adjust the P values by false-discovery rate
    merged$combined_fdr = p.adjust(merged$p_combined, method="BH", n=num_tests)
    
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
    clust = read.table(clust, header=TRUE, sep="\t")
    enrich = read.table(enrich, header=TRUE, sep="\t")
    meta_clust = read.table(meta_clust, header=TRUE, sep="\t")
    meta_enrich = read.table(meta_enrich, header=TRUE, sep="\t")
    
    meta = combine_enrichment_and_clustering(meta_enrich, meta_clust)
    ddd = combine_enrichment_and_clustering(enrich, clust)
    
    # read in p values from HPO similarity analyses
    phenotypes = read.table(phenotype, header=TRUE, sep="\t")
    
    ddd = merge(ddd, phenotypes, by="hgnc", all.x=TRUE)
    ddd$p_min_with_pheno = apply(ddd[, c("p_combined", "hpo_similarity_p_value")], 1, fishersMethod)
    
    # need names that are more informative as same across files, add prefix
    names(ddd)[3:length(names(ddd))] = paste("ddd", names(ddd)[3:length(names(ddd))], sep=".")
    names(meta)[3:length(names(meta))] = paste("meta", names(meta)[3:length(names(meta))], sep=".")
    
    # merge together files, focusing on genes with DNMs in DDD
    merged = merge(meta, ddd, by=c("hgnc", "chrom"), all.x=TRUE)
    merged$p_min = apply(merged[, c("ddd.p_combined", "ddd.p_min_with_pheno", "meta.p_combined")], 1,
        function(x) ifelse(length(x[!is.na(x)]) > 0, min(x[!is.na(x)]), NA))
    
    return(merged)
}
