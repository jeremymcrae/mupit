# function to open de novo data from the DDD study
#

#' get standardised de novo data for DDD study.
#'
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
standardise_ddd_de_novos <- function() {
    
    # load a set of DDD de novos, that are the independent event (ie duplicate
    # de novos within families have been removed).
    variants = read.table(file.path("data-raw", "de_novo_datasets",
        "de_novos.ddd_4k.ddd_only.txt"), header=TRUE,
        sep="\t", stringsAsFactors=FALSE, comment.char="")
    
    # standardise the SNV or INDEL flag
    variants$type = "indel"
    variants$type[variants$var_type == "DENOVO-SNP"] = "snv"
    
    # standardise the columns, and column names
    variants$person_id = gsub(" ", "", variants$person_stable_id)
    variants$chrom = gsub(" ", "", variants$chrom)
    variants$start_pos = gsub(" ", "", variants$pos)
    variants$ref_allele = gsub(" ", "", variants$ref)
    variants$alt_allele = gsub(" ", "", variants$alt)
    variants$end_pos = as.character(as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1)
    
    variants$hgnc = variants$symbol
    
    variants$study_code = "ddd_unpublished"
    variants$publication_doi = NA
    variants$study_phenotype = "developmental_disorders"
    variants$sex = variants$gender
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
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
get_ddd_de_novos <- function(diagnosed=NULL, subset=NULL) {
    
    # try and use the data from the package, if it is available, otherwise
    # generate the dataset
    variants = standardise_ddd_de_novos()
    
    # remove diagnosed patients, if maximising power
    if (!is.null(diagnosed)) {
        variants = variants[!(variants$person_id %in% diagnosed$id), ]
    }
    
    # sometimes we only want to use a subset of the DDD, such as when we are
    # investigating a single disease type, e.g. seizures. Then we will have
    # specified the DDD sample IDs that we wish to restrict to.
    if (!is.null(subset)) {
        variants = variants[variants$person_id %in% subset, ]
    }
    
    return(variants)
}
