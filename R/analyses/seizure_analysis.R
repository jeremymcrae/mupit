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
    
    samples = read.table(phenotype_filename, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # find the samples with either "seizure" or "epilep" in their HPO terms. 
    # this probably should check the HPO tree below the top term of "seizures", 
    # but this should suffice for now.
    seizure = grepl("[Ss]eizures", samples$child_terms)
    epilepsy = grepl("[Ep]pilep", samples$child_terms)
    has_seizure_term = seizure | epilepsy
    samples = samples[has_seizure_term, ]
    
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
#' @param meta true/false for whether to include meta-analysis populations
#' 
#' @return data frame containing HGNC, chrom, position, consequence, SNV or INDEL
#'     type, and study ID.
get_de_novos <- function(diagnosed, meta=FALSE) {
    
    seizure_probands = get_ddd_probands_with_seizures(PHENOTYPE_PATH, ID_MAP_PATH) 
    sample_ids = seizure_probands$proband_id
    
    variants = get_ddd_de_novos(diagnosed, subset=sample_ids)
    
    if (meta) {
        # read in other datasets and calculate numbers of LoF and NS, SNVs and indels
        seizure_de_novos = published_de_novos[published_de_novos$study_phenotype == "epilepsy", ]
        variants = rbind(variants, seizure_de_novos)
    }
    
    return(variants)
}

#' defines the cohort sizes, used to get the overall population size
#' 
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' 
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(diagnosed, meta=FALSE) {
    
    seizure_probands = get_ddd_probands_with_seizures(PHENOTYPE_PATH, ID_MAP_PATH)
    
    # number of trios studied in our data, minus the samples with diagnoses, in 
    # order to maximise the power
    undiagnosed = !(seizure_probands$proband_id %in% diagnosed$id)
    male = sum(seizure_probands$sex[undiagnosed] == "Male") # male probands
    female = sum(seizure_probands$sex[undiagnosed] == "Female") # female probands
    
    if (meta) {
        # sum up males and females across studies
        male = male + sum(cohorts$unique_male[cohorts$study_phenotype == "epilepsy"])
        female = female + sum(cohorts$unique_female[cohorts$study_phenotype == "epilepsy"])
    }
    
    return(list(male=male, female=female))
}

main <- function() {
    # diagnosed_path = "/nfs/users/nfs_j/jm33/apps/mupit/data-raw/Diagnoses_1133.txt"
    # diagnosed = get_ddd_diagnosed(diagnosed_path)
    diagnosed = get_likely_diagnosed()
    
    # analyse the DDD only de novos
    trios = get_trio_counts(diagnosed)
    de_novos = get_de_novos(diagnosed)
    enriched = analyse_gene_enrichment(de_novos, trios, "results/de_novos.seizures.ddd_only.manhattan.pdf")
    
    # analyse the DDD + other cohorts de novos
    trios_meta = get_trio_counts(diagnosed, meta=TRUE)
    de_novos_meta = get_de_novos(diagnosed, meta=TRUE)
    enriched_meta = analyse_gene_enrichment(de_novos_meta, trios_meta, "results/de_novos.seizures.meta-analysis.manhattan.pdf")
    
    write.table(enriched, file=file.path("results", 
        "de_novos.seizures.ddd_only.enrichment_results.txt"), sep="\t", 
        row.names=FALSE, quote=FALSE)
    
    write.table(enriched_meta, file=file.path("results", 
        "de_novos.seizures.meta-analysis.enrichment_results.txt"), sep="\t", 
        row.names=FALSE, quote=FALSE)
    
    head(enriched[order(enriched$p_func), ])
}


main()


