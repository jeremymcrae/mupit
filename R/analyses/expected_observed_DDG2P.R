# R code to estimate expected and observed numbers of mutations in DDG2P genes
# for different mutation categories.

library(mupit)

DDG2P_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/DDG2P_freeze_with_gencode19_genomic_coordinates_20141118_fixed.txt"

#' get the observed mutation counts in each of the gene sets
#'
#' Currently the counts are hard coded, it might be a good idea to convert
#' this to loading them dynamically, so that it can be more easily adapted to
#' different datasets
#'
#' @return list of lists defining the mutation counts in each gene set, under
#'     each mutation category.
get_observed_values <- function(de_novos, ddg2p, rates) {
    
    # restrict the de novos to ones ion genes for which we have mutation rate
    # data available, otherwise the missing genes will give an excess
    de_novos = de_novos[de_novos$hgnc %in% rates$hgnc, ]
    
    # find which de novos are in autosomal recessive genes, based on their
    # presence in the DDG2P
    index = get_ddg2p_index(ddg2p, de_novos)
    
    # convert the VEP consequences to one of five categories
    category = rep(NA, nrow(de_novos))
    category[de_novos$consequence %in% c("stop_gained", "splice_acceptor_variant",
        "splice_donor_variant")] = "lof"
    category[de_novos$consequence %in% c("missense_variant",
    "initiator_codon_variant", "stop_lost", "splice_region_variant")] = "missense"
    category[de_novos$consequence %in% c("synonymous_variant")] = "silent"
    category[de_novos$consequence %in% c("inframe_deletion", "inframe_insertion")] = "inframe"
    category[de_novos$consequence %in% c("frameshift_variant")] = "frameshift"
    
    de_novos$category = category
    
    # count each of the mutation types for the different gene lists
    ddg2p_ar = as.list(table(de_novos$category[index$ar]))
    ddg2p_non_ar = as.list(table(de_novos$category[index$non_ar]))
    nonddg2p = as.list(table(de_novos$category[!index$in_ddg2p]))
    
    return(list(ddg2p_ar=ddg2p_ar, ddg2p_non_ar=ddg2p_non_ar, nonddg2p=nonddg2p))
}

#' gets index positions for gene subsets of the rates datasets
#'
#' @param dataset dataframe, containing "hgnc" column, so that we can match HGNC
#'     symbols.
#' @param rates_data the Daly mutation rates data frame
#'
#' @return a dataframe of indexes
get_ddg2p_index <- function(ddg2p, dataset) {
    
    # define the inheritance modes for autosomal recessive and non-autosomal
    # recessive
    autosomal_recessive = c("Biallelic")
    non_autosomal_recessive_chrX = c("Hemizygous", "X-linked dominant")
    non_autosomal_recessive = c("Both", "Monoallelic")
    non_autosomal_recessive = c(non_autosomal_recessive, non_autosomal_recessive_chrX)
    
    # figure out different subsets of the DDG2P dataset
    # get gene names for ar and non.ar inheritance modes
    autosomal_recessive_genes = ddg2p$ddg2p_gene_name[ddg2p$Allelic_requirement %in% autosomal_recessive]
    non_autosomal_recessive_genes = ddg2p$ddg2p_gene_name[ddg2p$Allelic_requirement %in% non_autosomal_recessive]
    
    # find the genes from the rates data which are autosomal recessive etc
    autosomal_recessive = dataset$hgnc %in% autosomal_recessive_genes
    non_autosomal_recessive = dataset$hgnc %in% non_autosomal_recessive_genes
    in_ddg2p = autosomal_recessive | non_autosomal_recessive
    
    # join the indexes
    indexes = data.frame(hgnc=dataset$hgnc, non_ar=non_autosomal_recessive,
        ar=autosomal_recessive, in_ddg2p=in_ddg2p)
    
    return(indexes)
}

#' calculate the expected number of mutations for each gene
#'
#' @param rates list containing vectors of mutation rates for each gene under
#'     different mutation categories.
#' @param index vector of positions indicating the subset of genes to use (for
#'     example to indicate the DDG2P genes, or the autosomal recessive set).
#' @param inverse: whether to use the inverse of the index positions.
#'
#' @return dataframe of expected numbers of mutations for different mutation
#'     categories.
calculate_expected_numbers <- function(rates, index, inverse=FALSE) {
    
    # if we want the inverse, swap the index positions to refer to the
    # opposite positions
    if (inverse) { index = !index }
    
    # the expected values are the sums of rates for each mutation category
    lof = sum(rates$snv.lof.rate[index])
    missense = sum(rates$snv.missense.rate[index])
    silent = sum(rates$snv.silent.rate[index])
    inframe = sum(rates$indel.missense.rate[index])
    frameshift = sum(rates$indel.lof.rate[index])
    
    expected = data.frame(lof=lof, missense=missense, silent=silent,
        inframe=inframe, frameshift=frameshift)
    
    return(expected)
}

#' find the chance of the observed mutation count, given the expected values
#'
#' Uses a Poisson model to estimate the chance of finding as many mutations
#' as were observed, given the expected number of mutations (determined from
#' the mutation rates).
#'
#' @param observed data frame of counts of observed mutations, for different mutation
#'     mutation categories
#' @param expected data frame  of counts of expected mutations, for different
#'     mutation categories
#'
#' @return data frame of P values for each of the mutation categories.
calculate_p_from_observed <- function(observed, expected) {
    
    # NOTE: I'm not sure why the observed count is decremented by one.
    lof = ppois(observed$lof - 1, expected$lof, FALSE)
    missense = ppois(observed$missense - 1, expected$missense, FALSE)
    silent = ppois(observed$silent - 1, expected$silent, FALSE)
    inframe = ppois(observed$inframe - 1, expected$inframe, FALSE)
    frameshift = ppois(observed$frameshift - 1, expected$frameshift, FALSE)
    
    p_values = data.frame(lof=lof, missense=missense, silent=silent,
        inframe=inframe, frameshift=frameshift)
    
    return(p_values)
}

main <- function() {
    
    # read in DDG2P, for only Confirmed and Probable genes
    ddg2p = read.table(DDG2P_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    ddg2p = ddg2p[ddg2p$DDD_category != "Possible DD Gene", ]
    
    # get the mutation rates
    rates = get_mutation_rates(trios)
    
    # number of trios studied in our data
    trios = list()
    trios$male = 2408  # trios with male offspring
    trios$female = 1887  # trios with female offspring
    
    observed = get_observed_values(ddd_de_novos, ddg2p, rates)
    
    # find which genes have autosomal recessive (or not) inheritance
    indexes = get_ddg2p_index(ddg2p, rates)
    
    # calculate the expected number of mutated genes (for each mutation
    # category) for each set of genes
    expected_ddg2p_ar = calculate_expected_numbers(rates, indexes$ar)
    expected_ddg2p_non_ar = calculate_expected_numbers(rates, indexes$non_ar)
    expected_nonddg2p = calculate_expected_numbers(rates, indexes$in_ddg2p, inverse=TRUE)
    
    # estimate P values for each of the gene sets, across each mutation category
    ddg2p_ar_p = calculate_p_from_observed(observed$ddg2p_ar, expected_ddg2p_ar)
    ddg2p_non_ar_p = calculate_p_from_observed(observed$ddg2p_non_ar, expected_ddg2p_non_ar)
    nonddg2p_p = calculate_p_from_observed(observed$nonddg2p, expected_nonddg2p)
    
    # mimic the original output from this script
    t(expected_ddg2p_ar)
    t(expected_ddg2p_non_ar)
    t(expected_nonddg2p)
    
    t(ddg2p_ar_p)
    t(ddg2p_non_ar_p)
    t(nonddg2p_p)
}


main()
