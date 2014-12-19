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
    male_ddd = 2408 # male probands
    female_ddd = 1887 # female probands
    
    return(list(male = male_ddd, female = female_ddd))
}

main <- function() {
    # here's an example of how to use the functions in this script
    diagnosed_path = "/nfs/users/nfs_j/jm33/apps/mupit/data-raw/Diagnoses_1133.txt"
    diagnosed = get_ddd_diagnosed(diagnosed_path)
    num = get_trio_counts(diagnosed)
    num.trios.male = num$male
    num.trios.female = num$female
    
    # open the de novos, and 
    de_novos = get_ddd_de_novos(diagnosed)
    
    enriched = analyse_gene_enrichment(de_novos, num.trios.male, num.trios.female)
    head(enriched[order(enriched$p.func.daly), ])
}


main()



