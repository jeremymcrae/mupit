#### Calculate Daly mutation rates #######
# 
# Loads mutations rates from the Daly dataset, and adjusts for sex-chromosome 
# transmissions.

#' correct mutations rates for sex-chromosome transmission rates
#' 
#' @param rates list of the Daly dataframe, along with vectors of mutation
#'         rates for genes under different mutation classes
#' @param males number of trios with male offspring
#' @param females number of trios with female offspring
#' @export
#' 
#' @return list containing mutation rates, where chrX genes have been adjusted
#'     for the sex-specific transmission rates.
correct_for_x_chrom <- function(rates, males, females) {
    
    # TODO: why is male.transmissions euqal to the number of females? Why is it
    # TODO: not equal to the number of males?
    autosomal.transmissions = 2 * (males + females)
    female.transmissions = males + females
    male.transmissions = females
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish 
    # Family Health Study) phased de novo data.
    alpha = 3.4 
    male.chrx.scaling = 2 / (1 + (1 / alpha))
    female.chrx.scaling = 2 / (1 + alpha)
    
    # correct the non-PAR chrX genes for fewer transmissions and lower rate 
    # (dependent on alpha)
    chrX = which(rates$chrom == "X")
    x_scale_factor = ((male.transmissions * male.chrx.scaling) + (female.transmissions * female.chrx.scaling))/autosomal.transmissions
    
    rates$snv.silent.rate[chrX] = rates$snv.silent.rate[chrX] * x_scale_factor
    rates$snv.missense.rate[chrX] = rates$snv.missense.rate[chrX] * x_scale_factor
    rates$snv.lof.rate[chrX] = rates$snv.lof.rate[chrX] * x_scale_factor
    rates$indel.missense.rate[chrX] = rates$indel.missense.rate[chrX] * x_scale_factor
    rates$indel.lof.rate[chrX] = rates$indel.lof.rate[chrX] * x_scale_factor
    
    return(rates)
}

#' adapt indel rates for lower rate estimate from validated de novos
#' 
#' @param rates data frame of mutation rates.
#' @export
#' 
#' @return the rates data frame, with adjusted indel rates.
adjust_indel_rates <- function(rates) {
    
    # I think the following numbers were derived from the Daly dataset
    valid.nonsense = 102
    valid.frameshift = 95
    indel.scaling = 1.25 # ratio of frameshift to nonsense
    
    rates$indel.missense.rate = (rates$indel.missense.rate / indel.scaling) * (valid.frameshift / valid.nonsense)
    rates$indel.lof.rate = (rates$indel.lof.rate / indel.scaling) * (valid.frameshift / valid.nonsense)
    
    return(rates)
}

#' gets mutation rates from the Daly dataset
#' 
#' @param trios list of male and female proband counts in the dataset
#' @export
#' 
#' @return a list containing mutation rates for genes under different mutation 
#'     classes, along with the original Daly mutation rate dataset.
get_mutation_rates <- function(trios) {
    
    num.trios = trios$male + trios$female 
    auto.transmissions = 2 * num.trios
    
    daly = read.delim(file.path("data-raw", "fixed_mut_prob_fs_adjdepdiv.txt"), header=TRUE)
    
    # add chromosome annotation
    daly = merge(daly, gene_info, by.x="gene", by.y="hgnc", all.x=T) 
    
    # get the number of expected mutations, given the number of transmissions
    rates = data.frame(hgnc = daly$gene, chrom = daly$chrom)
    rates$snv.silent.rate = (10^daly$syn) * auto.transmissions
    rates$snv.missense.rate = (10^daly$mis + 10^daly$rdt) * auto.transmissions
    rates$snv.lof.rate = (10^daly$non + 10^daly$css) * auto.transmissions
    rates$indel.missense.rate = ((10^daly$frameshift) / 9) * auto.transmissions
    rates$indel.lof.rate = (10^daly$frameshift) * auto.transmissions

    # could scale to take account of longer transcript
    
    rates = adjust_indel_rates(rates)
    
    # catch cases where there is no rdt or css mutation rate, resulting in an 
    # NA for composite rates
    rates$snv.missense.rate[is.na(rates$snv.missense.rate)] = (10^daly$mis[is.na(rates$snv.missense.rate)]) * auto.transmissions
    rates$snv.lof.rate[is.na(rates$snv.lof.rate)] = (10^daly$non[is.na(rates$snv.lof.rate)]) * auto.transmissions
    
    # and correct for the X-chromosome rates
    rates = correct_for_x_chrom(rates, trios$male, trios$female)
    
    return(rates)
}

