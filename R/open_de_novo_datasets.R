# function to open de novo data from the DDD study
# 

#' find diagnosed probands in the DDD study, to exclude them from our data
#' 
#' @param path path to file defining diagnosed probands
#' 
#' @export
#' @return A list containing vectors with DECIPHER IDs, and sex of the diagnosed
#'     probands
get_ddd_diagnosed <- function(path) {
    
    # read in samples that have been diagnosed, so as to remove from our data
    diagnoses = read.delim(path, header=TRUE)
    index = rowSums(diagnoses[, c(4:14)]) > 0
    
    diagnosed = list()
    diagnosed$id = diagnoses$DDD_ID[index]
    diagnosed$sex = diagnoses$Sex[index]
    
    return(diagnosed)
}

#' get standardised de novo data for DDD study.
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
standardise_ddd_de_novos <- function() {
    
    variants = read.table(file.path("data-raw", "de_novo_datasets", 
        "DNG_Variants_20Feb2014_NonRed_Clean_NoTwins_NoClusters.txt"), header=TRUE, 
        sep="\t", stringsAsFactors=FALSE, comment.char="")
    
    # remove diagnosed patients, if maximising power
    variants = variants[-which(variants$DECIPHER_ID %in% diagnosed$id), ]
    
    # sometimes we only want to use a subset of the DDD, such as when we are 
    # investigating a single disease type, e.g. seizures. Then we will have
    # specified the DDD sample IDs that we wish to restrict to.
    if (!is.null(sample_ids)[1]) {
        variants = variants[variants$person_stable_id %in% sample_ids, ]
    }
    
    # standardise the SNV or INDEL flag
    variants$type = "indel"
    variants$type[variants$snp_or_indel == "DENOVO-SNP"] = "snv"
    
    # standardise the columns, and column names
    variants$person_id = gsub(" ", "", variants$person_stable_id)
    variants$chrom = gsub(" ", "", variants$chr)
    variants$start_pos = gsub(" ", "", variants$pos)
    variants$ref_allele = gsub(" ", "", variants$ref)
    variants$alt_allele = gsub(" ", "", variants$alt)
    variants$end_pos = as.character(as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1)
    
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "ddd_unpublished"
    variants$publication_doi = NA
    variants$study_phenotype = "developmental_disorders"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence", 
        "study_code", "publication_doi", "study_phenotype", "type"))
    
    return(variants)
}

#' get de novo data for DDD study.
#' 
#' @param diagnosed dataframe of samples in DDD with a diagnosis, who should be
#'     excluded from analyses
#' @param subset vector of sample IDs, if we want to restrict to a specific
#'         set of probands, such as when we are analysing a phenotype-specific
#'         subset of the DDD.
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
get_ddd_de_novos <- function(diagnosed, subset=NULL) {
    
    # try and use the data from the package, if it is available, otherwise 
    # generate the dataset
    if (exists("ddd_de_novos") & nrow(ddd_de_novos) > 0) {
        variants = ddd_de_novos
     } else {
        variants = standardise_ddd_de_novos()
        ddd_de_novos = variants
        save(ddd_de_novos, file="data/ddd_de_novos.rda", compress="xz")
    }
    
    # remove diagnosed patients, if maximising power
    variants = variants[!(variants$person_id %in% diagnosed$id), ]
    
    # sometimes we only want to use a subset of the DDD, such as when we are 
    # investigating a single disease type, e.g. seizures. Then we will have
    # specified the DDD sample IDs that we wish to restrict to.
    if (!is.null(subset)[1]) {
        variants = variants[variants$person_id %in% subset, ]
    }
    
    return(variants)
}
