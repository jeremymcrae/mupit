# script to analyse enrichment of de novos in genes in probands with seizures
# in the DDD

library(mupit)

#' defines the cohort sizes, used to get the overall population size
#' 
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' 
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(diagnosed) {
    
    # number of trios studied in our data
    male_ddd = 582 # trios with male offspring
    female_ddd = 548 # trios with female offspring
    
    # remove diagnosed patients, if maximising power
    male_ddd = male_ddd - length(which(diagnosed$sex == "Male"))
    female_ddd = female_ddd - length(which(diagnosed$sex == "Female"))
    
    male = sum(cohorts$unique_male, male_ddd)
    female = sum(cohorts$unique_female, female_ddd)
    
    return(list(male = male, female = female))
}

#' combine datasets listing de novo mutations into a single data frame
#' 
#' @param diagnosed list of IDs and sex for probands with diagnoses in the DDD
#' 
#' @return data frame containing HGNC, chrom, position, consequence, SNV or 
#'    INDEL type, and study ID.
open_datasets <- function(diagnosed) {
    
    ddd = open_ddd_de_novos(diagnosed)
    data = rbind(ddd, published_de_novos)
    
    return(data)
}

main <- function() {
    # here's an example of how to use the functions in this script
    diagnosed_path = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis/data-raw/Diagnoses_1133.txt"
    diagnosed = get_ddd_diagnosed(diagnosed_path)
    num = get_trio_counts(diagnosed)
    num.trios.male = num$male
    num.trios.female = num$female
    
    # open the de novos, and 
    de_novos = open_data sets(diagnosed)
    
    enriched = analyse_gene_enrichment(de_novos, num.trios.male, num.trios.female)
    head(enriched[order(enriched$p.func.daly), ])
}


main()



