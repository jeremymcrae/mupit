# Calculate mutation rates
#
# Loads gene-based mutation rates, in order to determine the expected number of
# mutations per gene, given the number of studied probands and adjusts for
# sex-chromosome transmissions.

#' gets mutation rates from the gene-based datasets
#'
#' This defaults to the gene-based mutation rates from Nature Genetics
#' 46:944-950 (2014) doi:10.1038/ng.3050, but we can pass in other gene-based
#' mutation rate datasets.
#'
#' @param trios list of male and female proband counts in the dataset
#' @param rates optionally specify a dataframe containing per-gene mutation
#'    rates, or NULL
#' @export
#'
#' @return a dataframe of mutation rates for genes under different mutation
#'     classes.
get_gene_based_mutation_rates <- function(trios, rates=NULL) {
    
    # Get the number of autosomal transmissions, so we can estimate the expected
    # number of mutations given the number of potential transmissions.
    autosomal = 2 * (trios$male + trios$female)
    
    # if we haven't passed in a rates dataset, default to the gene rates from
    # Samocha et al., Nature Genetics 46:944-950.
    if (is.null(rates)) { rates = mupit::gene_rates }
    
    # add chromosome annotation to the rates, if it isn't already included,
    # so that we can later correct for chrX transmission differences
    if (!("chrom" %in% names(rates))) {
        rates = merge(rates, mupit::gene_info, by="hgnc", all.x=TRUE)
        rates$chrom = as.character(rates$chrom)
        rates$chrom[is.na(rates$chrom)] = "UNKNOWN"
    }
    
    # convert rates from log-scaled values, so we can multiply by the number of
    # transmissions
    columns = c("syn", "mis", "non", "splice_site", "frameshift")
    rates[, columns] = 10 ** rates[, columns]
    
    # Some genes lack splice site mutation rates (they have NA instead). This
    # impedes estimating the snv.lof.rate, which is estimated from the nonsense
    # and splice_site rates. Convert these NA splice_site values to zero.
    rates$splice_site[is.na(rates$splice_site)] = 0
    
    # get the number of expected mutations, given the number of transmissions
    rates$snv.silent.rate = rates$syn * autosomal
    rates$snv.missense.rate = rates$mis * autosomal
    rates$snv.lof.rate = (rates$non + rates$splice_site) * autosomal
    rates$indel.missense.rate = (rates$frameshift / 9) * autosomal
    rates$indel.lof.rate = rates$frameshift * autosomal
    
    # correct for the known ratio of indels to nonsense, and for transmissions
    # on the X-chromosome
    rates = adjust_indel_rates(rates)
    rates = correct_for_x_chrom(rates, trios$male, trios$female)
    
    # subset to the columns we need to estimate enrichment probabilities
    rates = rates[, c("hgnc", "chrom", "snv.silent.rate", "snv.missense.rate",
        "snv.lof.rate", "indel.missense.rate", "indel.lof.rate")]
    
    return(rates)
}

#' calculate mutation rates based on the gene length
#'
#' This has been discontinued in favor of using rates derived from de novo
#' mutation rates, from Nature Genetics 46:944-950.
#'
#' @param trios list of male and female proband counts in the dataset
#' @export
#'
#' @return a dataframe of mutation rates for genes under different mutation
#'     classes.
get_length_based_rates <- function(trios) {
    
    cds.length = mupit::gene_info$cds_length
    
    # Get the number of autosomal transmissions, so we can estimate the expected
    # number of mutations given the number of potential transmissions.
    autosomal = 2 * (trios$male + trios$female)
    
    # set mutation rates for snvs and indels
    snv_rate = 1.5E-8 # higher than genome-wide mutation rate, due to higher GC ...
    # TODO: check that the following indel mutation rate is correct, since the
    # TODO: comment claims it should be 10% of the SNV mutation rate, but it
    # TODO: currently looks  like 30%
    indel_rate = 0.53E-9 # ~10% of genome-wide SNV mutation rate, no reason to think higher in exome
    
    # specify proportion of coding mutations of different types
    props = list()
    props$snv.lof = 0.0485 # from Daly
    props$snv.missense = 0.6597 # from Daly
    props$snv.silent = 0.2918 # from Daly
    
    props$indel.lof = 0.9 # non-3n, from size distribution in neutral sequence
    props$indel.missense = 0.1
    
    # calculate rates of missense and lof mutations, multiply by 2 for two
    # transmissions per child and number of trios
    rates = data.frame(hgnc = mupit::gene_info$hgnc, chrom = mupit::gene_info$chrom)
    rates$snv.silent.rate = cds.length * snv_rate * props$snv.silent * autosomal
    rates$snv.missense.rate = cds.length * snv_rate * props$snv.missense * autosomal
    rates$snv.lof.rate = cds.length * snv_rate * props$snv.lof * autosomal
    rates$indel.missense.rate = cds.length * indel_rate * props$indel.missense * autosomal
    rates$indel.lof.rate = cds.length * indel_rate * props$indel.lof * autosomal
    
    rates = correct_for_x_chrom(rates, trios$male, trios$female)
    
    return(rates)
}

#' correct mutations rates for sex-chromosome transmission rates
#'
#' @param rates gene-based data frame, containing rates for different mutation
#'     classes.
#' @param male_n number of trios with male offspring
#' @param female_n number of trios with female offspring
#' @export
#'
#' @return a dataframe of mutation rates for genes under different mutation
#'     classes.
correct_for_x_chrom <- function(rates, male_n, female_n) {
    
    # figure out the number of transmissions for autosomal, male and female
    # transmissions
    # TODO: why is male.transmissions equal to the number of females? Why is it
    # TODO: not equal to the number of males?
    autosomal = 2 * (male_n + female_n)
    female = male_n + female_n
    male = female_n
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    alpha = 3.4
    male_factor = 2 / (1 + (1 / alpha))
    female_factor = 2 / (1 + alpha)
    
    # correct the non-PAR chrX genes for fewer transmissions and lower rate
    # (dependent on alpha)
    chrX = rates$chrom == "X"
    x_factor = ((male * male_factor) + (female * female_factor))/autosomal
    
    rates$snv.silent.rate[chrX] = rates$snv.silent.rate[chrX] * x_factor
    rates$snv.missense.rate[chrX] = rates$snv.missense.rate[chrX] * x_factor
    rates$snv.lof.rate[chrX] = rates$snv.lof.rate[chrX] * x_factor
    rates$indel.missense.rate[chrX] = rates$indel.missense.rate[chrX] * x_factor
    rates$indel.lof.rate[chrX] = rates$indel.lof.rate[chrX] * x_factor
    
    return(rates)
}

#' adapt indel rates for lower rate estimate from validated de novos
#'
#' The indel mutation rates from Samocha et al., Nature Genetics 46:944-950
#' assume that the overall indel mutation rate is 1.25-fold greater than the
#' overall nonsense mutation rate, ie there are 1.25 times as many frameshifts
#' as nonsense mutations. We have our own estimates for the ratio, derived from
#' our de novo validation efforts, which we shall apply in place of the Samocha
#' et al ratios.
#'
#' @param rates data frame of mutation rates.
#' @export
#'
#' @return the rates data frame, with adjusted indel rates.
adjust_indel_rates <- function(rates) {
    
    # I think the following numbers were derived from the DDD dataset
    valid_nonsense = 102
    valid_frameshift = 95
    ddd_ratio = valid_frameshift / valid_nonsense
    samocha_ratio = 1.25  # Nature Genetics 46:944-950 frameshift to nonsense ratio
    
    rates$indel.missense.rate = (rates$indel.missense.rate / samocha_ratio) * ddd_ratio
    rates$indel.lof.rate = (rates$indel.lof.rate / samocha_ratio) * ddd_ratio
    
    return(rates)
}
