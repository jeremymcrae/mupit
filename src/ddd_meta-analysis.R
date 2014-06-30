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

get_trio_counts <- function(diagnosed) {
    # defines the cohort sizes, used to get the overall population size
    # 
    # Args:
    #     diagnosed: list of sex and ID for probands diagnosed in the DDD study
    # 
    # Returns:
    #     list with total counts of trios with male and female offspring
    
    # number of trios studied in our data
    male.ddd = 582 # trios with male offspring
    female.ddd = 548 # trios with female offspring
    
    # remove diagnosed patients, if maximising power
    male.ddd = male.ddd - length(which(diagnosed$sex == "Male"))
    female.ddd = female.ddd - length(which(diagnosed$sex == "Female"))
    
    # number of trios studied in deligt data
    male.deligt = 47
    female.deligt = 53
    
    # number of trios studied in autism data
    male.autism = 764
    female.autism = 183
    
    # number of trios studied in rauch data
    male.rauch = 32
    female.rauch = 19
    
    # number of trios studied in fromer data
    male.fromer = 317
    female.fromer = 306
    
    # number of trios studied in epi4k data
    male.epi4k = 156
    female.epi4k = 108
    
    # number of trios studied in zaidi data
    male.zaidi = 220
    female.zaidi = 142
    
    # sum up males and females across studies
    male = male.ddd + male.deligt + male.autism + male.rauch + male.fromer + male.epi4k + male.zaidi
    female = female.ddd + female.deligt + female.autism + female.rauch + female.fromer + female.epi4k + female.zaidi
    
    return(list(male = male, female = female))
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
    
    ddd = open_ddd_denovos(diagnosed)
    
    # read in other datasets and calculate numbers of LoF and NS, SNVs and indels
    rauch = open_rauch_de_novos()
    deligt = open_deligt_de_novos()
    epi4k = open_epi4k_de_novos()
    autism = open_autism_de_novos()
    fromer = open_fromer_de_novos ()
    zaidi = open_zaidi_de_novos ()
    
    data = rbind(ddd, rauch, deligt, epi4k, autism, fromer, zaidi)
    
    return(data)
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



