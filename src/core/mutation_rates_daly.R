#### Calculate Daly mutation rates #######
# 
# Loads mutations rates from the Daly dataset, and adjusts for sex-chromosome 
# transmissions.

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")

correct_for_x_chrom <- function(rates, gene_info, males, females) {
    # correct mutations rates for sex-chromosome transmission rates
    # 
    # Args:
    #     rates: list of the Daly dataframe, along with vectors of mutation
    #         rates for genes under different mutation classes
    #     males: number of trios with male offspring
    #     females: number of trios with female offspring
    # 
    # Returns:
    #     list containing mutation rates, where chrX genes have been adjusted
    #     for the sex-specific transmission rates.
    
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
    chrX = which(gene_info$chr == "X")
    x_scale_factor = ((male.transmissions * male.chrx.scaling) + (female.transmissions * female.chrx.scaling))/autosomal.transmissions
    
    rates$snv.silent.rate[chrX] = rates$snv.silent.rate[chrX] * x_scale_factor
    rates$snv.missense.rate[chrX] = rates$snv.missense.rate[chrX] * x_scale_factor
    rates$snv.lof.rate[chrX] = rates$snv.lof.rate[chrX] * x_scale_factor
    rates$indel.missense.rate[chrX] = rates$indel.missense.rate[chrX] * x_scale_factor
    rates$indel.lof.rate[chrX] = rates$indel.lof.rate[chrX] * x_scale_factor
    
    return(rates)
}

adjust_indel_rates <- function(rates) {
    # adapt indel rates for lower rate estimate from validated de novos
    # 
    # Args:
    #     rates: data frame of mutation rates.
    # 
    # Returns:
    #     the rates data frame, with adjusted indel rates.
    
    # I think the following numbers were derived from the Daly dataset
    valid.nonsense = 102
    valid.frameshift = 95
    indel.scaling = 1.25 # ratio of frameshift to nonsense
    
    rates$indel.missense.rate = (rates$indel.missense.rate / indel.scaling) * (valid.frameshift / valid.nonsense)
    rates$indel.lof.rate = (rates$indel.lof.rate / indel.scaling) * (valid.frameshift / valid.nonsense)
    
    return(rates)
}

get_mutation_rates <- function(num.trios.male, num.trios.female) {
    # gets mutation rates from the Daly dataset
    # 
    # Args:
    #     num.trios.male: number of trios with a male offspring
    #     num.trios.female: number of trios with a female offspring
    # 
    # Returns:
    #     a list containing mutation rates for genes under different mutation 
    #     classes, along with the original Daly mutation rate dataset.
    
    num.trios = num.trios.male + num.trios.female 
    auto.transmissions = 2 * num.trios
    
    gene.info = read.delim(file.path(DATA_DIR, "CDS_LENGTH_B37_chr.txt"), header=T)
    daly = read.delim(file.path(DATA_DIR, "fixed_mut_prob_fs_adjdepdiv.txt"), header=T)
    
    # add chromosome annotation
    daly = merge(daly, gene.info, by.x=2, by.y=4, all.x=T) 
    
    # get the number of expected mutations, given the number of transmissions
    rates = data.frame(HGNC = daly$gene)
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
    rates = correct_for_x_chrom(rates, daly, num.trios.male, num.trios.female)
    
    values = list(rates = rates, gene.info = daly)
    
    return(values)
}

