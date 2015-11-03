# program to calculate the significance of seeing N DNMs of a specific
# combination of functional types in a particular gene in M trios

# RATIONALE: use gene coding sequence to predict rate of DNMs in coding sequence
# for each gene of different functional classes, then estimate the probability
# of seeing the observed combination of different functional classes of DNMs
# assuming number of DNMs in each class is Poisson distributed.


#' tallies the mutation types observed for each gene
#'
#' @param de_novos data frame listing all the de novo mutations, with columns
#'        for HGNC symbol, consequence type (VEP style predictions), and a
#'        column indicating SNV, or indel.
#' @export
#'
#' @return data frame with tally of de novo mutations for each of the mutation
#'         types
#'
#' @examples
#' vars = read.table(header=TRUE, text="
#'    person_id   hgnc     chrom   start_pos   consequence        type
#'    person_1    ARID1B   6       157431695   missense_variant   snv
#'    person_2    ARID1B   6       157502190   stop_gained        snv")
#' get_de_novo_counts(vars)
get_de_novo_counts <- function(de_novos) {
    
    de_novos$hgnc = as.character(de_novos$hgnc)
    de_novos$consequence = as.character(de_novos$consequence)
    
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant", "initiator_codon_variant", "start_lost",
        "conserved_exon_terminus_variant")
    missense_cq = c("missense_variant", "stop_lost", "inframe_deletion",
        "inframe_insertion", "coding_sequence_variant", "protein_altering_variant")
    
    # group the lof and missence consequence strings, and drop all the
    # non-functional de novos
    de_novos$consequence[de_novos$consequence %in% lof_cq] = "lof"
    de_novos$consequence[de_novos$consequence %in% missense_cq] = "missense"
    de_novos = de_novos[de_novos$consequence %in% c("missense", "lof"), ]
    
    # count the number of de novos for each type/consequence combination
    de_novo_counts = reshape::cast(de_novos, hgnc ~ consequence + type, value = "person_id", length)
    
    # include the positions of the de novos at the minimum position for each gene
    by_hgnc = split(de_novos[, c("hgnc", "chrom", "start_pos")], de_novos$hgnc)
    de_novo_counts$chrom = sapply(by_hgnc, function(x) x[["chrom"]][1])
    de_novo_counts$min_pos = sapply(by_hgnc, function(x) (sort(as.numeric(x[["start_pos"]])))[1])
    
    # ensure all the required mutation categories are available as columns
    for (column in c("lof_indel", "missense_indel", "lof_snv", "missense_snv")) {
        if (!column %in% names(de_novo_counts)) { de_novo_counts[[column]] = 0 }
    }
    
    return(data.frame(de_novo_counts))
}

#' tests whether genes are enriched with de novo mutations
#'
#' @param rates gene mutation rates per consequence type
#' @param counts data frame with tally of de novo mutations per gene for each of
#'     the mutation types: lof_snv, lof_indel, missense_snv, missense_indel.
#' @param all_genes whether to test all genes in the genome (most will test
#'     the probability of observing 0 de novos).
#'
#' @export
#' @return data frame with gene info, mutation rates and P values from testing
#'     for enrichment.
#'
#' @examples
#' counts = read.table(header=TRUE, text="
#'     hgnc   lof_indel lof_snv missense_indel missense_snv chrom min_pos
#'     ARID1B 1       1            1         1              6     157150547
#'     KMT2A  0       0            0         1              11    118367083")
#'
#' rates = read.table(header=TRUE, text="
#'     hgnc   chrom snv.lof.rate indel.lof.rate snv.missense.rate indel.missense.rate
#'     ARID1B 6     0.01        0.005          0.05             0.005
#'     KMT2A  11    0.01        0.005          0.05             0.005")
#'
#' gene_enrichment(rates, counts)
gene_enrichment <- function(rates, counts, all_genes=FALSE) {
    
    observed = merge(counts, rates, by = c("hgnc", "chrom"), all.x=TRUE)
    
    # occasionally we want results for all genes, not just the ones we have
    # observed de novos at (basically just testing the probably observing 0
    # de novos in the absent genes). We use this full set for plotting QQ plots
    if (all_genes) {
        observed = merge(counts, rates, by = c("hgnc", "chrom"), all=TRUE)
        observed[is.na(observed)] = 0
    }
    
    # for each gene, sum the de novo counts across SNVs and indels for the
    # different functional categories: synonymous, loss of function, missense
    # and functional
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
    observed$p_lof = ppois(lof_count - 1, lambda=lof_rate, lower.tail=FALSE)
    observed$p_func = ppois(func_count - 1, lambda=func_rate, lower.tail=FALSE)
    
    # remove the rate columns from the dataframe
    observed = observed[, !grepl("rate", names(observed))]
    
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
#'
#' @examples
#' trios = list(male=1000, female=1000)
#' vars = read.table(header=TRUE, text="
#'      person_id hgnc chrom start_pos consequence type
#'      person_1 ARID1B 6 157150547 inframe_deletion indel
#'      person_2 ARID1B 6 157431695 missense_variant snv
#'      person_3 ARID1B 6 157454186 frameshift_variant indel
#'      person_4 ARID1B 6 157502190 stop_gained snv
#'      person_4 KMT2A 11 118367083 missense_variant snv")
#' rates = read.table(header=TRUE, text="
#'      hgnc    chrom  syn  mis  non   splice_site  frameshift
#'      ARID1B  6      -5   -6   -6.5  -5           -7
#'      KMT2A   11     -5   -6   -6.5  -5           -7")
#' analyse_gene_enrichment(vars, trios, rates=rates)
analyse_gene_enrichment <- function(de_novos, trios, plot_path=NULL, all_genes=FALSE, rates=NULL) {
    
    # tally the de novos by consequence and variant type
    de_novo_counts = get_de_novo_counts(de_novos)
    
    # get the mutation rates for each gene
    mutation_rates = get_gene_based_mutation_rates(trios, rates)
    
    # calculate p values for each gene using the mutation rates
    enriched = gene_enrichment(mutation_rates, de_novo_counts, all_genes)
    
    # make a manhattan plot of enrichment P values
    if (!is.null(plot_path)) {
        num_tests = 18500
        plot_enrichment_graphs(enriched, num_tests, plot_path)
    }
    
    # remove the position column (which is only used to be able to locate the
    # gene's position on a chromosome on a Manhattan plot).
    enriched$min_pos = NULL
    
    return(enriched)
}
