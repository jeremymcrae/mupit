# function to open de novo data from the DDD study
# 

#' get de novo data for DDD study.
#' 
#' @param diagnosed dataframe of samples in DDD with a diagnosis, who should be
#'     excluded from analyses
#' @param sample_ids vector of sample IDs, if we want to restrict to a specific
#'         set of probands, such as when we are analysing a phenotype-specific
#'         subset of the DDD.
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
open_ddd_de_novos <- function(diagnosed, sample_ids=NULL) {
    
    de_novos = read.delim(file.path("data-raw", "de_novo_datasets", "DNG_Variants_20Feb2014_NonRed_Clean_NoTwins_NoClusters.txt"), header=TRUE, colClasses = "character")
    
    # remove diagnosed patients, if maximising power
    de_novos = de_novos[-which(de_novos$DECIPHER_ID %in% diagnosed$id), ]
    
    # sometimes we only want to use a subset of the DDD, such as when we are 
    # investigating a single disease type, e.g. seizures. Then we will have
    # specified the DDD sample IDs that we wish to restrict to.
    if (!is.null(sample_ids)[1]) {
        de_novos = de_novos[de_novos$person_stable_id %in% sample_ids, ]
    }
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(de_novos$snp_or_indel == "DENOVO-SNP")
    de_novos$TYPE = "INDEL"
    de_novos$TYPE[TYPE.index] = "SNV"
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("curated_HGNC", "curated_CQ", "pos", "chr", "TYPE"))
    names(de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    de_novos$STUDY = "DDD"
    
    return(de_novos)
}
