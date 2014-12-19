# script to analyse enrichment of de novos in genes in probands with seizures
# in the DDD

library(mupit)

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"
PHENOTYPE_PATH = file.path(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")
ID_MAP_PATH = file.path(DATAFREEZE_DIR, "person_sanger_decipher.txt")

#' find which DDD probands have seizures, and return their IDs and sexes
#' 
#' @param phenotype_filename path to file containing phenotype information,
#'         including whether they have seizures
#' 
#' @return a data frame containing DDD sample IDs, their decipher IDs, and their
#'     sex
get_ddd_probands_with_seizures <- function(phenotype_filename, id_path) {
    
    samples = read.table(phenotype_filename, sep="\t", header=TRUE, colClasses="character")
    samples = samples[grepl("[Ss]eizures", samples$child_terms), ]
    
    ddd_ids = read.table(id_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    ddd_ids = ddd_ids[!grepl(":", ddd_ids$decipher_id), ]
    ddd_ids = ddd_ids[, c("person_stable_id", "decipher_id")]
    ddd_ids = ddd_ids[!duplicated(ddd_ids), ]
    
    samples = merge(samples, ddd_ids, by.x="patient_id", by.y="decipher_id")
    
    probands = data.frame(proband_id = samples$person_stable_id, decipher_id = 
        samples$patient_id, sex = samples$gender)
    
    return(probands)
}

#' combine datasets listing de novo mutations into a single data frame
#' 
#' @param diagnosed list of IDs and sex for probands with diagnoses in the DDD
#' 
#' @return data frame containing HGNC, chrom, position, consequence, SNV or INDEL
#'     type, and study ID.
get_de_novos <- function(diagnosed) {
    
    seizure_probands = get_ddd_probands_with_seizures(PHENOTYPE_PATH, ID_MAP_PATH) 
    sample_ids = seizure_probands$proband_id
    
    ddd_de_novos = get_ddd_de_novos(diagnosed, subset=sample_ids)
    
    # read in other datasets and calculate numbers of LoF and NS, SNVs and indels
    seizure_de_novos = published_de_novos[published_de_novos$study_phenotype == "epilepsy", ]
    variants = rbind(ddd_de_novos, seizure_de_novos)
    
    return(variants)
}

#' defines the cohort sizes, used to get the overall population size
#' 
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' 
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(diagnosed) {
    
    seizure_probands = get_ddd_probands_with_seizures(PHENOTYPE_PATH, ID_MAP_PATH)
    
    # number of trios studied in our data, minus the samples with diagnoses, in 
    # order to maximise the power
    undiagnosed = !(seizure_probands$proband_id %in% diagnosed$id)
    male_ddd = sum(seizure_probands$sex[undiagnosed] == "Male") # male probands
    female_ddd = sum(seizure_probands$sex[undiagnosed] == "Female") # female probands
    
    # sum up males and females across studies
    male = male_ddd + sum(cohorts$unique_male[cohorts$study_phenotype == "epilepsy"])
    female = female_ddd + sum(cohorts$unique_female[cohorts$study_phenotype == "epilepsy"])
    
    return(list(male = male, female = female))
}


main <- function() {
    diagnosed_path = "/nfs/users/nfs_j/jm33/apps/mupit/data-raw/Diagnoses_1133.txt"
    diagnosed = get_ddd_diagnosed(diagnosed_path)
    num = get_trio_counts(diagnosed)
    num.trios.male = num$male
    num.trios.female = num$female
    
    # open the de novos, and 
    de_novos = get_de_novos(diagnosed)
    
    enriched = analyse_gene_enrichment(de_novos, num.trios.male, num.trios.female)
    
    head(enriched[order(enriched$p.func.daly), ])
}


main()


