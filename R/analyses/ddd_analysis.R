# script to analyse enrichment of de novo mutations in genes in probands

library(mupit)

#' defines the cohort sizes, used to get the overall population size
#' 
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' @param meta true/false for whether to include meta-analysis populations
#' 
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(diagnosed, meta=FALSE) {
    
    # number of trios studied in our data
    male = 2408 # male probands
    female = 1887 # female probands
    
    # remove diagnosed patients, if maximising power
    male = male - sum(diagnosed$sex == "Male")
    female = female - sum(diagnosed$sex == "Female")
    
    if (meta) {
        male = male + sum(cohorts$unique_male)
        female = female + sum(cohorts$unique_female)
    }
    
    return(list(male=male, female=female))
}

#' combine datasets listing de novo mutations into a single data frame
#' 
#' @param diagnosed list of IDs and sex for probands with diagnoses in the DDD
#' @param meta true/false for whether to include meta-analysis populations
#' 
#' @return data frame containing HGNC, chrom, position, consequence, SNV or 
#'    INDEL type, and study ID.
get_de_novos <- function(diagnosed, meta=FALSE) {
    
    variants = get_ddd_de_novos(diagnosed)
    
    if (meta) {
        variants = rbind(variants, published_de_novos)
    }
    
    return(variants)
}

main <- function() {
    # here's an example of how to use the functions in this script
    diagnosed_path = "/nfs/users/nfs_j/jm33/apps/mupit/data-raw/Diagnoses_1133.txt"
    diagnosed = get_ddd_diagnosed(diagnosed_path)
    
    # analyse the DDD only de novos
    trios = get_trio_counts(diagnosed)
    de_novos = get_de_novos(diagnosed)
    enriched = analyse_gene_enrichment(de_novos, trios)
    
    # analyse the DDD+ other cohorts de novos
    trios_meta = get_trio_counts(diagnosed, meta=TRUE)
    de_novos_meta = get_de_novos(diagnosed, meta=TRUE)
    enriched_meta = analyse_gene_enrichment(de_novos_meta, trios_meta)
    
    write.table(enriched, file=file.path("results", 
        "de_novos.ddd.ddd_only.enrichment_results.txt"), sep="\t", 
        row.names=FALSE, quote=FALSE)
    
    write.table(enriched_meta, file=file.path("results", 
        "de_novos.ddd.meta-analysis.enrichment_results.txt"), sep="\t", 
        row.names=FALSE, quote=FALSE)
    
    head(enriched[order(enriched$p_func), ])
}


main()



