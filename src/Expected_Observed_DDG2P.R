# R code to estimate expected and observed numbers of mutations in DDG2P genes
# for different mutation categories.

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")
SRC_DIR = file.path(CODE_DIR, "src")
source(file.path(SRC_DIR, "mutation_rates_daly.R"))

get_observed_values <- function() {
    # get the observed mutation counts in each of the gene sets 
    # 
    # Currently the counts are hard coded, it might be a good idea to convert
    # this to loading them dynamically, so that it can be more easily adapted to
    # different datasets
    # 
    # Returns:
    #     list of lists defining the mutation counts in each gene set, under 
    #     each mutation category.
    
    # define the counts for autosomal recessive genes in the DDG2P
    ddg2p.ar = list()
    ddg2p.ar$lof = 6
    ddg2p.ar$missense = 43
    ddg2p.ar$silent = 12
    ddg2p.ar$inframe = 0
    ddg2p.ar$frameshift = 0
    
    # define the counts for non-autosomal recessive genes in the DDG2P
    ddg2p.non.ar = list()
    ddg2p.non.ar$lof = 45
    ddg2p.non.ar$missense = 97
    ddg2p.non.ar$silent = 9
    ddg2p.non.ar$inframe = 1
    ddg2p.non.ar$frameshift = 38
    
    # define the counts for non-DDG2P genes
    nonddg2p = list()
    nonddg2p$lof = 86
    nonddg2p$missense = 724
    nonddg2p$silent = 239
    nonddg2p$inframe = 9
    nonddg2p$frameshift = 57
    
    return(list(ddg2p.ar = ddg2p.ar, ddg2p.non.ar = ddg2p.non.ar, nonddg2p = nonddg2p))
}

get_ddg2p_index <- function(rates_data) {
    # gets index positions for gene subsets of the rates datasets
    # 
    # Args:
    #     rates_data: the Daly mutation rates data frame
    # 
    # Returns:
    #     a list of index position vectors
    
    # read in DDG2P, for only Confirmed and Probable genes
    ddg2p = read.delim(file.path(DATA_DIR, "DDG2P_Confirmed_Probable_20131107.txt"), header=T)
    
    # define the inheritance modes for autosomal recessive and non-autosomal 
    # recessive
    ar.CQ = c("Biallelic")
    non.ar.chrx.CQ = c("Hemizygous", "X-linked dominant")
    non.ar.auto.CQ = c("Both", "Monoallelic")
    non.ar.CQ = c(non.ar.auto.CQ, non.ar.chrx.CQ)
    
    # figure out different subsets of the DDG2P dataset
    # get gene names for ar and non.ar inheritance modes
    ddg2p.ar.uniq = unique(ddg2p$gene[ddg2p$mode %in% ar.CQ])
    ddg2p.non.ar.uniq = unique(ddg2p$gene[ddg2p$mode %in% non.ar.CQ]) 
    
    # make the indexes of positions
    ar.index = which(rates_data$gene %in% ddg2p.ar.uniq)
    non.ar.index = which(rates_data$gene %in% ddg2p.non.ar.uniq)
    ddg2p.index = unique(c(non.ar.index, ar.index))
    
    # join the indexes into a list
    indexes = list(non.ar.index = non.ar.index, ar.index = ar.index, 
        ddg2p.index = ddg2p.index)
    
    return(indexes)
}

calculate_expected_numbers <- function(rates, index, inverse=FALSE) {
    # calculate the expected number of mutations for each gene
    # 
    # Args:
    #     rates: list containing vectors of mutation rates for each gene under 
    #         different mutation categories.
    #     index: vector of positions indicating the subset of genes to use (for
    #         example to indicate the DDG2P genes, or the autosomal recessive 
    #         set).
    #     inverse: whether to use the inverse of the index positions.
    # 
    # Returns:
    #     list of expected numbers of mutations for different mutation 
    #     categories.
    
    # if we want the inverse, swap the index positions to refer to the 
    # opposite positions
    if (inverse) {
        index = !(seq(rates$snv.lof.rate) %in% index)
    }
    
    # the expected values are the sums of rates for each mutation category
    expected = list()
    expected$lof = sum(rates$snv.lof.rate[index])
    expected$missense = sum(rates$snv.missense.rate[index])
    expected$silent = sum(rates$snv.silent.rate[index])
    expected$inframe = sum(rates$indel.missense.rate[index])
    expected$frameshift = sum(rates$indel.lof.rate[index])
    
    return(expected)
}

calculate_p_from_observed <- function(observed, expected) {
    # find the chance of the observed mutation count, given the expected values
    # 
    # Uses a Poisson model to estimate the chance of finding as many mutations
    # as were observed, given the expected number of mutations (determined from
    # the mutation rates).
    # 
    # Args:
    #     observed: list of counts of observed mutations, for different mutation
    #         categories
    #     expected: list of counts of expected mutations, for different mutation
    #         categories
    # 
    # Returns:
    #     list of P values for each of the mutation categories.
    
    p_values = list()
    # NOTE: I'm not sure why the observed count is decremented by one.
    p_values$lof = ppois(observed$lof-1, expected$lof, F)
    p_values$missense = ppois(observed$missense-1, expected$missense, F)
    p_values$silent = ppois(observed$silent-1, expected$silent, F)
    p_values$inframe = ppois(observed$inframe-1, expected$inframe, F)
    p_values$frameshift = ppois(observed$frameshift-1, expected$frameshift, F)
    
    return(p_values)
}

main <- function(){
    
    observed = get_observed_values()
    
    # number of trios studied in our data
    num.trios.male = 582 # trios with male offspring
    num.trios.female = 548 # trios with female offspring
    
    # get the mutation rates
    rates = get_mutation_rates(num.trios.male, num.trios.female)
    daly = rates$daly
    
    # find which genes have autosomal recessive (or not) inheritance
    indexes = get_ddg2p_index(daly)
    ar.index = indexes$ar.index
    non.ar.index = indexes$non.ar.index
    ddg2p.index = indexes$ddg2p.index
    
    # calculate the expected number of mutated genes (for each mutation 
    # category) for each set of genes
    expected.ddg2p.ar = calculate_expected_numbers(rates, ar.index)
    expected.ddg2p.non.ar = calculate_expected_numbers(rates, non.ar.index)
    expected.nonddg2p = calculate_expected_numbers(rates, ddg2p.index, inverse=TRUE)
    
    # estimate P values for each of the gene sets, across each mutation category
    ddg2p.ar.p = calculate_p_from_observed(observed$ddg2p.ar, expected.ddg2p.ar)
    ddg2p.non.ar.p = calculate_p_from_observed(observed$ddg2p.non.ar, expected.ddg2p.non.ar)
    nonddg2p.p = calculate_p_from_observed(observed$nonddg2p, expected.nonddg2p)
    
    # mimic the original output from this script
    t(rbind(expected.ddg2p.ar))
    t(rbind(expected.ddg2p.non.ar))
    t(rbind(expected.nonddg2p))
    
    t(rbind(ddg2p.ar.p))
    t(rbind(ddg2p.non.ar.p))
    t(rbind(nonddg2p.p))
}


main()
