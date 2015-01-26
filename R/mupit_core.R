# program to calculate the significance of seeing N DNMs of a specific
# combination of functional types in a particular gene in M trios

# RATIONALE: use gene coding sequence to predict rate of DNMs in coding sequence
# for each gene of different functional classes, then estimate the probability
# of seeing the observed combination of different functional classes of DNMs
# assuming number of DNMs in each class is Poisson distributed.


#' tallies the mutation types observed for each gene
#'
#' @param de_novos data frame listing all the de novo mutations, with columns
#'     for HGNC symbol, consequence type (VEP style predictions), and a column
#'     indicating SNV, or indel.
#' @export
#'
#' @return data frame with tally of de novo mutations for each of the mutation
#'     types
get_de_novo_counts <- function(de_novos) {
    
    de_novos$hgnc = as.character(de_novos$hgnc)
    de_novos$consequence = as.character(de_novos$consequence)
    
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant")
    missense_cq = c("missense_variant", "initiator_codon_variant", "stop_lost",
        "inframe_deletion", "inframe_insertion", "splice_region_variant")
    synonymous = "synonymous_variant"
    
    lof_regex = paste(lof_cq, collapse = "|")
    missense_regex = paste(missense_cq, collapse = "|")
    
    # group the lof and missence consequence strings, and drop all the
    # non-functional de novos
    de_novos$consequence[grepl(lof_regex, de_novos$consequence)] = "lof"
    de_novos$consequence[grepl(missense_regex, de_novos$consequence)] = "missense"
    de_novos$consequence[grepl(synonymous, de_novos$consequence)] = "synonymous"
    de_novos = de_novos[grepl("missense|lof|synonymous", de_novos$consequence), ]
    
    # count the number of de novos for each type/consequence combination
    de_novo_counts = reshape::cast(de_novos, hgnc ~ consequence + type, value = "person_id", length)
    
    # include the positions of the de novos at the minimum position for each gene
    de_novos$min_pos = "min_pos"
    de_novos$temp = "chrom"
    pos = reshape::cast(de_novos, hgnc ~ min_pos, value = "start_pos", min)
    chrom = reshape::cast(de_novos, hgnc ~ temp, value = "chrom", min)
    de_novo_counts = merge(pos, de_novo_counts, by = "hgnc")
    de_novo_counts = merge(chrom, de_novo_counts, by = "hgnc")
    
    # ensure all the required mutation categories are available as columns
    if (!("lof_indel" %in% names(de_novo_counts))) { de_novo_counts$lof_indel = 0 }
    if (!("missense_indel" %in% names(de_novo_counts))) { de_novo_counts$missense_indel = 0 }
    if (!("lof_snv" %in% names(de_novo_counts))) { de_novo_counts$lof_snv = 0 }
    if (!("missense_snv" %in% names(de_novo_counts))) { de_novo_counts$missense_snv = 0 }
    if (!("synonymous_snv" %in% names(de_novo_counts))) { de_novo_counts$synonymous_snv = 0 }
    
    return(de_novo_counts)
}

#' tests whether genes are enriched with de novo mutations
#'
#' @param rates gene mutation rates per consequence type
#' @param counts data frame with tally of de novo mutations per gene for each of
#'     the mutation types: lof_snv, lof_indel, missense_snv, missense_indel.
#' @param num_tests number of tests performed (used for multiple correction).
#' @param all_genes whether to test all genes in the genome (most will test
#'     the probability of observing 0 de novos).
#'
#' @export
#' @return data frame with gene info, mutation rates and P values from testing
#'     for enrichment.
test_enrichment <- function(rates, counts, num_tests, all_genes=FALSE) {
    
    observed = merge(counts, rates, by = c("hgnc", "chrom"), all.x=TRUE)
    
    # occasionally we want results for all genes, not just the ones we have
    # observed de novos at (basically just testing the probably observing 0
    # de novos in the absent genes). We use this full set for plotting QQ plots,
    # since
    if (all_genes) {
        observed = merge(counts, rates, by = c("hgnc", "chrom"), all=TRUE)
        observed[is.na(observed)] = 0
        
        if (num_tests < nrow(observed)) { num_tests = nrow(observed) }
    }
    
    # for each gene, sum the de novo counts across SNVs and indels for the
    # different functional categories: synonymous, loss of function, missense
    # and functional
    synonymous_count = observed$synonymous_snv
    lof_count = observed$lof_snv + observed$lof_indel
    missense_count = observed$missense_snv + observed$missense_indel
    func_count = lof_count + missense_count
    
    # for each gene, sum the mutation rates across SNVs and indels for the
    # different functional categories
    synonymous_rate = observed$snv.silent.rate
    lof_rate = observed$snv.lof.rate + observed$indel.lof.rate
    missense_rate = observed$snv.missense.rate + observed$indel.missense.rate
    func_rate = lof_rate + missense_rate
    
    # calculate the probably of getting the observed number of de novos, given
    # the mutation rate
    observed$p_synonymous = dpois(synonymous_count, lambda=synonymous_rate)
    observed$p_lof = dpois(lof_count, lambda=lof_rate)
    observed$p_func = dpois(func_count, lambda=func_rate)
    
    # correct the P values for multiple testing by false discovery rate
    observed$fdr_synonymous = p.adjust(observed$p_lof, method="BH", n=num_tests)
    observed$fdr_lof = p.adjust(observed$p_lof, method="BH", n=num_tests)
    observed$fdr_func = p.adjust(observed$p_func, method="BH", n=num_tests)
    
    return(observed)
}

#' analyse whether de novo mutations are enriched in genes
#'
#' @param de_novos data frame containing all the observed de novos for all the
#'     genes
#' @param trios list of male and female proband counts in the population
#' @param plot_path path to save enrichment plots to
#' @param all_genes whether to test all genes in the genome (most will test
#'     the probability of observing 0 de novos).
#' @param rates gene-based mutation rates data frame, or NA
#' @export
#'
#' @return data frame containing results from testiong for enrichment of de
#'     in each gene with de novos in it.
analyse_gene_enrichment <- function(de_novos, trios, plot_path=NA, all_genes=FALSE, rates=NA) {
    
    # tally the de novos by consequence and variant type
    de_novo_counts = get_de_novo_counts(de_novos)
    
    # get the mutation rates for each gene
    mutation_rates = get_gene_based_mutation_rates(trios, rates)
    
    # calculate p values for each gene using the mutation rates
    num_tests = 18500
    enriched = test_enrichment(mutation_rates, de_novo_counts, num_tests, all_genes)
    
    # make a manhattan plot of enrichment P values
    if (!is.na(plot_path)) {
        plot_enrichment_graphs(enriched, num_tests, plot_path)
    }
    
    return(enriched)
}
