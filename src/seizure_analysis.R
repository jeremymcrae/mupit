# script to analyse enrichment of de novos in genes in probands with seizures
# in the DDD


CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")
SRC_DIR = file.path(CODE_DIR, "src")
DE_NOVO_DIR = file.path(DATA_DIR, "de_novo_datasets")
source(file.path(SRC_DIR, "core", "mup-it_undiagnosed.R"))

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/"
PHENOYTYPE_FILENAME = file.path(DATAFREEZE_DIR, "phenotypes.shared.pcs.relatedness.diagnosis.20140415.txt")

CQ.LOF = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant")
CQ.MISSENSE = c("missense_variant", "initiator_codon_variant", "stop_lost", "inframe_deletion", "inframe_insertion")

get_ddd_probands_with_seizures <- function(phenotype_filename) {
    
    samples = read.table(phenotype_filename, sep = "\t", header = TRUE, colClasses = "character")
    samples = samples[samples$seizures == 1, ]
    
    decipher_id = samples$decipher_id
    sample_id = samples$DDD_ID
    sex = samples$Sex
    
    probands = list(sample_id = sample_id, sex = sex, decipher_id = decipher_id)
    
    return(probands)
}

open_datasets <- function(diagnosed) {
    # combine datasets listing de novo mutations into a single data frame
    # 
    # Args:
    #     diagnosed: list of IDs and sex for probands with diagnoses in the DDD
    # 
    # Returns:
    #     data frame containing HGNC, chrom, position, consequence, SNV or INDEL
    #     type, and study ID.
    
    seizure_probands = get_ddd_probands_with_seizures(PHENOYTYPE_FILENAME) 
    sample_ids = seizure_probands$sample_id
    
    ddd = open_ddd_denovos(diagnosed, sample_ids)
    
    # read in other datasets and calculate numbers of LoF and NS, SNVs and indels
    epi4k = open_epi4k_de_novos()
    
    data = rbind(ddd, epi4k)
    
    return(data)
}

get_trio_counts <- function(diagnosed) {
    # defines the cohort sizes, used to get the overall population size
    # 
    # Args:
    #     diagnosed: list of sex and ID for probands diagnosed in the DDD study
    # 
    # Returns:
    #     list with total counts of trios with male and female offspring
    
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
    diagnosed = get_ddd_diagnosed()
    num = get_trio_counts(diagnosed)
    num.trios.male = num$male
    num.trios.female = num$female
    
    # open the de novos, and 
    de_novos = open_datasets(diagnosed)
    
    analyse_gene_enrichment(de_novos, num.trios.male, num.trios.female)
}


main()


