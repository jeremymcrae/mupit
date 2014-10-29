# Calculate mutation rates based on the length of transcripts for a gene
# 
# This has been discontinued in favor of using rates derived from de novo 
# mutation rates, as provided by Mark Daly


#' calculate mutation rates based on the gene length
#' 
#' @param num.trios.male number of trios with male offspring in the dataset
#' @param num.trios.female number of trios with female offspring in the dataset
#' @export
#' 
#' @return list containing vectors of mutation rates
get_length_based_rates <- function(num.trios.male, num.trios.female) {
    
    cds.length = gene_info$cds_length
    chrX = which(gene_info$chrom == "X")
    
    auto.transmissions = 2 * (num.trios.male + num.trios.female)
    
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
    
    props$indel.lof = 0.9 # non-3n, from size distribution in neutral sequence
    props$indel.missense = 0.1
    
    # calculate rates of missense and lof mutations, multiply by 2 for two 
    # transmissions per child and number of trios
    rates = data.frame(hgnc = gene_info$hgnc)
    rates$snv.missense.rate = cds.length * snv_rate * props$snv.missense * auto.transmissions
    rates$snv.lof.rate = cds.length * snv_rate * props$snv.lof * auto.transmissions
    rates$indel.missense.rate = cds.length * indel_rate * props$indel.missense * auto.transmissions
    rates$indel.lof.rate = cds.length * indel_rate * props$indel.lof * auto.transmissions
    
    rates = correct_length_rate_for_x_chrom(rates, num.trios.male, num.trios.female, cds.length, chrX, snv_rate, indel_rate, props)
    
    return(rates)
}

#' correct mutations rates for sex-chromosome transmission rates
#' 
#' @param rates dataframe of mutation rates for different mutation categories
#' @param males number of trios with male offspring
#' @param females number of trios with female offspring
#' @param cds.length vector of cds lengths for genes, in same order as rates 
#'         data frame
#' @param chrX index of row numbers for rates data frame, where the gene is on 
#'         the X chromosome
#' @param snv_rate mutation rate for SNVs
#' @param indel_rate mutation rate for indels
#' @param props list of proportions of mutations as LOF SNVs, missense SNVs, 
#'         LOF indels and missense indels
#' @export
#' 
#' @return list containing mutation rates, where chrX genes have been adjusted
#'     for the sex-specific transmission rates.
correct_length_rate_for_x_chrom <- function(rates, males, females, cds.length, chrX, snv_rate, indel_rate, props) {
    
    # TODO: why is male.transmissions equal to the number of females? Why is it
    # TODO: not equal to the number of males?
    female.transmissions = males + females
    male.transmissions = females
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish 
    # Family Health Study) phased de novo data.
    alpha = 3.4 
    male.chrx.scaling = 2 / (1 + (1 / alpha))
    female.chrx.scaling = 2 / (1 + alpha)
    
    # set sex-specific mutation rates
    male.snv_rate = snv_rate * male.chrx.scaling
    female.snv_rate = snv_rate * female.chrx.scaling
    
    male.indel_rate = indel_rate * male.chrx.scaling
    female.indel_rate = indel_rate * female.chrx.scaling
    
    # correct non-PAR chrX genes for fewer transmissions and lower rate 
    # (dependent on alpha) currently doing for all chrX genes, not just non-PAR 
    # genes
    x_snv_scaling = (male.transmissions * male.snv_rate + female.transmissions * female.snv_rate)
    x_indel_scaling = (male.transmissions * male.indel_rate + female.transmissions * female.indel_rate)
    
    rates$snv.missense.rate[chrX] = cds.length[chrX] * props$snv.missense * x_snv_scaling
    rates$snv.lof.rate[chrX] = cds.length[chrX] * props$snv.lof * x_snv_scaling
    rates$indel.missense.rate[chrX] = cds.length[chrX] * props$indel.missense * x_indel_scaling
    rates$indel.lof.rate[chrX] = cds.length[chrX] * props$indel.lof * x_indel_scaling
    
    return(rates)
}
