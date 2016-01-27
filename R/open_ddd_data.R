# function to open de novo data from the DDD study
#

#' get standardised de novo data for DDD study.
#'
#' @param path path to DDD de novo dataset
#'
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
standardise_ddd_de_novos <- function(path) {
    
    variants = read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # standardise the SNV or INDEL flag
    variants$type = "indel"
    variants$type[variants$var_type == "DENOVO-SNP"] = "snv"
    
    # standardise the columns, and column names
    variants$person_id = variants$person_stable_id
    variants$start_pos = variants$pos
    variants$ref_allele = variants$ref
    variants$alt_allele = variants$alt
    variants$end_pos = variants$start_pos + nchar(variants$ref_allele) - 1
    
    variants$hgnc = variants$symbol
    
    variants$study_code = "ddd_unpublished"
    variants$publication_doi = NA
    variants$study_phenotype = "developmental_disorders"
    
    # standardise the sex codes
    variants$sex[variants$sex == "M"] = "male"
    variants$sex[variants$sex == "F"] = "female"
    
    variants = variants[, c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype", "type")]
    
    return(variants)
}

#' get de novo data for DDD study.
#'
#' @param path path to DDD de novo dataset
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
get_ddd_de_novos <- function(path, diagnosed=NULL, subset=NULL) {
    
    variants = standardise_ddd_de_novos(path)
    
    # remove diagnosed patients, if maximising power
    if (!is.null(diagnosed)) {
        variants = variants[!(variants$person_id %in% diagnosed$person_id), ]
    }
    
    # sometimes we only want to use a subset of the DDD, such as when we are
    # investigating a single disease type, e.g. seizures. Then we will have
    # specified the DDD sample IDs that we wish to restrict to.
    if (!is.null(subset)) {
        variants = variants[variants$person_id %in% subset, ]
    }
    
    return(variants)
}

#' load the current DDG2P dataset.
#'
#' @param path path to DDD de novo dataset
#'
#' @export
#' @return data frame of DDG2P genes
load_ddg2p <- function(path) {
    # the current file is missing a field from the header
    ddg2p = read.table(path, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
    
    # fix the header problem
    cols = ddg2p[1, ]
    cols[8] = "description"
    names(ddg2p) = as.character(cols)
    ddg2p = ddg2p[-1, ]
    
    # restrict outrselves to the high-confidence genes with a dominant mode of
    # inheritance
    ddg2p = ddg2p[ddg2p$type != "Possible DD Gene", ]
    ddg2p$dominant = ddg2p$mode %in% c("Monoallelic", "X-linked dominant")
    ddg2p$hemizygous = ddg2p$mode == "Hemizygous"
    
    # make sure the position columns are as numbers and clean the row numbers
    # after removing less confident genes
    ddg2p$start = as.numeric(ddg2p$start)
    ddg2p$stop = as.numeric(ddg2p$stop)
    row.names(ddg2p) = 1:nrow(ddg2p)
    
    return(ddg2p)
}
