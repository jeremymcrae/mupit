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
    
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant")
    missense_cq = c("missense_variant", "initiator_codon_variant", "stop_lost",
        "inframe_deletion", "inframe_insertion")
    
    lof_regex = paste(lof_cq, collapse = "|")
    missense_regex = paste(missense_cq, collapse = "|")
    all_regex = paste(c(lof_cq, missense_cq), collapse = "|")
    
    # count the number of de novos for each type/consequence combination
    counts = reshape::cast(de_novos, hgnc ~ type + consequence, value = "study", length)
    indels = counts[, grep("INDEL", names(counts))]
    snvs = counts[, grep("SNV", names(counts))]
    
    # and sum the de novo counts in each type/consequence category
    de_novo_counts = data.frame(hgnc = counts$hgnc)
    de_novo_counts$lof.snvs = rowSums(data.frame(snvs[, grep(lof_regex, names(snvs))]))
    de_novo_counts$missense.snvs = rowSums(data.frame(snvs[, grep(missense_regex, names(snvs))]))
    de_novo_counts$lof.indels = rowSums(data.frame(indels[, grep(lof_regex, names(indels))]))
    de_novo_counts$missense.indels = rowSums(data.frame(indels[, grep(missense_regex, names(indels))]))
    
    return(de_novo_counts)
}

#' tests whether genes are enriched with de novo mutations
#' 
#' @param rates gene mutation rates per consequence type
#' @param de_novos data frame containing all the observed de novos for all the 
#'            genes
#' @param counts data frame with tally of de novo mutations for each of the
#'            mutation types.
#' @param num.tests number of tests performed (used for multiple correction).
#' @export
#' 
#' @return data frame with gene info, mutation rates and P values from testing
#'     for enrichment.
get_p_values <- function(rates, de_novos, counts, num.tests) {
    
    # set-up vectors to store gene-specific information only for observed genes
    observed = data.frame(matrix(NA, nrow = nrow(counts), ncol = 8))
    names(observed) = c("chr", "coord", "snv.missense.rate", 
        "snv.lof.rate", "indel.missense.rate", "indel.lof.rate", "p.func", 
        "p.lof")

    # loop for each observed gene, test for functional variants and lof 
    # variants, using the gene mutation rates
    for (i in 1:nrow(counts)) {
        gene = as.character(counts$hgnc[i])
        
        # continue to next gene if mutation rates not available for the gene
        if (!(gene %in% rates$hgnc)) { next }
        
        gene.index = which(rates$hgnc == gene)
        
        # get the mutation rates for the gene
        observed[i, c("snv.missense.rate", "snv.lof.rate", "indel.missense.rate", "indel.lof.rate")] = rates[gene.index, c("snv.missense.rate", "snv.lof.rate", "indel.missense.rate", "indel.lof.rate")]
        
        # figure out the chromosome and nucleotide position of the gene
        data.index = which(de_novos$hgnc == gene)[1]
        observed$chr[i] = de_novos$chrom[data.index]
        observed$coord[i] = de_novos$position[data.index]
        
        # count the observed de novos in each functional category
        lof_count = sum(counts[i, c("lof.snvs", "lof.indels")])
        missense_count = sum(counts[i, c("missense.snvs", "missense.indels")])
        func_count = lof_count + missense_count
        
        # get the mutation rates for each functional category
        lof_rate = sum(observed[i, c("snv.lof.rate", "indel.lof.rate")])
        missense_rate = sum(observed[i, c("snv.missense.rate", "indel.missense.rate")]) 
        func_rate = lof_rate + missense_rate
        
        # calculate the probability of observing said de novos, given the 
        # gene mutation rates
        observed$p.lof[i] = dpois(lof_count, lambda=lof_rate)
        observed$p.func[i] = dpois(func_count, lambda=func_rate)
        
        if (i %% 100 == 0) {
            print(paste(i, " out of ", nrow(counts), " genes", sep = ""))
        }
    }
    
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
    p_vals_length = get_p_values(cds_rates, de_novos, de_novo_counts, num.tests)
    p_vals_daly = get_p_values(daly_rates, de_novos, de_novo_counts, num.tests)
    
    # write out results table
    enriched = cbind(de_novo_counts, p_vals_length, p_vals_daly)
    names(enriched) = c("hgnc", "LOF.snvs", "NS.snvs", "LOF.indels",
        "NS.indels", "chr", "coord", "snv.missense.rate",
        "snv.lof.rate", "indel.missense.rate", "indel.lof.rate", "p.func",
        "p.lof", "fdr.lof", "fdr.func", "chr", "coord",
        "daly.snv.missense.rate", "daly.snv.lof.rate",
        "daly.indel.missense.rate", "daly.indel.lof.rate", "daly.p.func",
        "daly.p.lof", "daly.fdr.lof", "daly.fdr.func")
    
    plot_enrichment_graphs(enriched, num.tests)
    
    return(enriched)
}



