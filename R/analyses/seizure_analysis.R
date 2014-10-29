# script to analyse enrichment of de novos in genes in probands with seizures
# in the DDD

library(mupit)

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/"
PHENOYTYPE_FILENAME = file.path(DATAFREEZE_DIR, "phenotypes.shared.pcs.relatedness.diagnosis.20140415.txt")

#' find which DDD probands have seizures, and return their IDs and sexes
#' 
#' @param phenotype_filename path to file containing phenotype information,
#'         including whether they have seizures
#' 
#' @return a data frame containing DDD sample IDs, their decipher IDs, and their
#'     sex
get_ddd_probands_with_seizures <- function(phenotype_filename) {
    
    samples = read.table(phenotype_filename, sep = "\t", header = TRUE, colClasses = "character")
    samples = samples[samples$seizures == 1, ]
    
    probands = data.frame(sample_id = samples$DDD_ID, decipher_id = 
        samples$decipher_id, sex = samples$Sex)
    
    return(probands)
}

#' combine datasets listing de novo mutations into a single data frame
#' 
#' @param diagnosed list of IDs and sex for probands with diagnoses in the DDD
#' 
#' @return data frame containing HGNC, chrom, position, consequence, SNV or INDEL
#'     type, and study ID.
open_datasets <- function(diagnosed) {
    
    seizure_probands = get_ddd_probands_with_seizures(PHENOYTYPE_FILENAME) 
    sample_ids = seizure_probands$sample_id
    
    ddd = open_ddd_de_novos(diagnosed, sample_ids)
    
    # read in other datasets and calculate numbers of LoF and NS, SNVs and indels
    data = rbind(ddd, epi4k_de_novos)
    
    return(data)
}

#' defines the cohort sizes, used to get the overall population size
#' 
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' 
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(diagnosed) {
    
    seizure_probands = get_ddd_probands_with_seizures(PHENOYTYPE_FILENAME)
    
    # number of trios studied in our data
    male.ddd = table(seizure_probands$sex)[["Male"]] # trios with male offspring
    female.ddd = table(seizure_probands$sex)[["Female"]] # trios with female offspring
    
    # remove diagnosed patients, if maximising power
    diagnosed_index = which(seizure_probands$decipher_id %in% diagnosed$id)
    male.ddd = male.ddd - length(which(seizure_probands$sex[diagnosed_index] == "Male"))
    female.ddd = female.ddd - length(which(seizure_probands$sex[diagnosed_index] == "Female"))
    
    # number of trios studied in epi4k data
    male.epi4k = 156
    female.epi4k = 108
    
    # sum up males and females across studies
    male = male.ddd + male.epi4k 
    female = female.ddd + female.epi4k
    
    return(list(male = male, female = female))
}


main <- function() {
    diagnosed_path = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis/data-raw/Diagnoses_1133.txt"
    diagnosed = get_ddd_diagnosed(diagnosed_path)
    num = get_trio_counts(diagnosed)
    num.trios.male = num$male
    num.trios.female = num$female
    
    # open the de novos, and 
    de_novos = open_datasets(diagnosed)
    
    enriched = analyse_gene_enrichment(de_novos, num.trios.male, num.trios.female)
    
    head(enriched[order(enriched$p.func.daly), ])
}


main()


