# program to calculate the significance of seeing N DNMs of a specific
# combination of functional types in a particular gene in M trios

# RATIONALE: use gene coding sequence to predict rate of DNMs in coding sequence
# for each gene of different functional classes, then estimate the probability
# of seeing the observed combination of different functional classes of DNMs
# assuming number of DNMs in each class is Poisson distributed.


#' tallies the mutation types observed for each gene
#' 
#' @param de_novos data frame listing all the de novo mutations, with columns
#'         for HGNC symbol, consequence type (VEP style predictions), and a
#'         column indicating SNV, or indel.
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
        "inframe_deletion", "inframe_insertion")
    
    lof_regex = paste(lof_cq, collapse = "|")
    missense_regex = paste(missense_cq, collapse = "|")
    
    # group the lof and missence consequence strings, and drop all the 
    # non-functional de novos
    de_novos$consequence[grepl(lof_regex, de_novos$consequence)] = "lof"
    de_novos$consequence[grepl(missense_regex, de_novos$consequence)] = "missense"
    de_novos = de_novos[grepl("missense|lof", de_novos$consequence), ]
    
    # count the number of de novos for each type/consequence combination
    de_novo_counts = reshape::cast(de_novos, hgnc ~ consequence + type, value = "person_id", length)
    
    # include the positions of the de novos at the minimum position for each gene
    de_novos$min_pos = "min_pos"
    de_novos$temp = "chrom"
    pos = reshape::cast(de_novos, hgnc ~ min_pos, value = "start_pos", min)
    chrom = reshape::cast(de_novos, hgnc ~ temp, value = "chrom", min)
    de_novo_counts = merge(pos, de_novo_counts, by = "hgnc")
    de_novo_counts = merge(chrom, de_novo_counts, by = "hgnc")
    
    return(de_novo_counts)
}

#' tests whether genes are enriched with de novo mutations
#' 
#' @param rates gene mutation rates per consequence type
#' @param counts data frame with tally of de novo mutations per gene for each of
#'     the mutation types: lof_snv, lof_indel, missense_snv, missense_indel.
#' @param num.tests number of tests performed (used for multiple correction).
#' 
#' @export
#' @return data frame with gene info, mutation rates and P values from testing
#'     for enrichment.
get_p_values <- function(rates, counts, num.tests) {
    
    observed = merge(counts, rates, by = c("hgnc", "chrom"), all.x=TRUE)
    
    # for each gene, sum the de novo counts across SNVs and indels for the
    # different functional categories: loss of function, missense and functional
    lof_count = observed$lof_snv + observed$lof_indel
    missense_count = observed$missense_snv + observed$missense_indel
    func_count = lof_count + missense_count
    
    # for each gene, sum the mutation rates across SNVs and indels for the
    # different functional categories
    lof_rate = observed$snv.lof.rate + observed$indel.lof.rate
    missense_rate = observed$snv.missense.rate + observed$indel.missense.rate
    func_rate = lof_rate + missense_rate
    
    # calculate the probably of getting the observed number of de novos, given
    # the mutation rate
    observed$p.lof = dpois(lof_count, lambda=lof_rate)
    observed$p.func = dpois(func_count, lambda=func_rate)
    
    # correct the P values for multiple testing by false discovery rate
    observed$fdr.lof = p.adjust(observed$p.lof, method="BH", n=num.tests)
    observed$fdr.func = p.adjust(observed$p.func, method="BH", n=num.tests)
    
    return(observed)
}

#' run the analysis of whether de novo mutations are enriched in genes
#' 
#' @param de_novos data frame containing all the observed de novos for all the 
#'         genes
#' @param num.trios.male number of trios with male offspring in the dataset
#' @param num.trios.female number of trios with female offspring in the dataset
#' @export
#' 
#' @return data frame containing results from testiong for enrichment of de
#'     in each gene with de novos in it.
analyse_gene_enrichment <- function(de_novos, num.trios.male, num.trios.female) {
    
    # tally the de novos by consequence and variant type
    de_novo_counts = get_de_novo_counts(de_novos)
    
    # get the length-based, and Daly-based mutation rates for each gene
    cds_rates = get_length_based_rates(num.trios.male, num.trios.female)
    daly_rates = get_mutation_rates(num.trios.male, num.trios.female)
    
    # calculate p values for each gene using the different mutation rates
    num.tests = 18500
    p_vals_length = get_p_values(cds_rates, de_novo_counts, num.tests)
    p_vals_daly = get_p_values(daly_rates, de_novo_counts, num.tests)
    
    # write out results table
    enriched = merge(p_vals_length, p_vals_daly, by = c("hgnc", "chrom", 
        "min_pos", "lof_indel", "lof_snv", "missense_indel", "missense_snv"))
    
    # fix the column names
    names(enriched) = c("hgnc", "chrom", "min_pos", "lof_indel", 
        "lof_snv", "missense_indel", "missense_snv", 
        "snv.missense.rate.length", "snv.lof.rate.length", 
        "indel.missense.rate.length", "indel.lof.rate.length", "p.lof.length", 
        "p.func.length", "fdr.lof.length", "fdr.func.length", 
        "snv.silent.rate.daly", "snv.missense.rate.daly", "snv.lof.rate.daly", 
        "indel.missense.rate.daly", "indel.lof.rate.daly", "p.lof.daly", 
        "p.func.daly", "fdr.lof.daly", "fdr.func.daly")
    
    plot_enrichment_graphs(enriched, num.tests)
    
    return(enriched)
}



