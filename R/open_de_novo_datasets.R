# functions to open de novo data from different datasets, and standardise the
# format, so that datasets can be combined
# 
# 

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")
SRC_DIR = file.path(CODE_DIR, "src")
DE_NOVO_DIR = file.path(DATA_DIR, "de_novo_datasets")

#' get de novo data for DDD study.
#' 
#' @param diagnosed
#' @param sample_ids vector of sample IDs, if we want to restrict to a specific
#'         set of probands, such as when we are analysing a phenotype-specific
#'         subset of the DDD.
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
open_ddd_denovos <- function(diagnosed, sample_ids=NA) {
    
    de_novos = read.delim(file.path(DE_NOVO_DIR, "DNG_Variants_20Feb2014_NonRed_Clean_NoTwins_NoClusters.txt"), header=TRUE, colClasses = "character")
    
    # remove diagnosed patients, if maximising power
    de_novos = de_novos[-which(de_novos$DECIPHER_ID %in% diagnosed$id), ]
    
    # sometimes we only want to use a subset of the DDD, such as when we are 
    # investigating a single disease type, e.g. seizures. Then we will have
    # specified the DDD sample IDs that we wish to restrict to.
    if (!is.na(sample_ids)[1]) {
        de_novos = de_novos[de_novos$person_stable_id %in% sample_ids, ]
    }
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(de_novos$snp_or_indel == "DENOVO-SNP")
    de_novos$TYPE = "INDEL"
    de_novos$TYPE[TYPE.index] = "SNV"
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("curated_HGNC", "curated_CQ", "pos", "chr", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "DDD"
    
    return(de_novos)
}

#' get de novo data for Rauch et al. intellectual disability exome study
#' 
#' De novo mutation data sourced from supplementary tables 2 and 3 from
#' Rauch et al. (2012) Lancet 380:1674–1682 
#' doi: 10.1016/S0140-6736(12)61480-9
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
open_rauch_de_novos <- function() {
   
    de_novos = read.delim(file.path(DE_NOVO_DIR, "rauch_v2.txt"), header=TRUE, colClasses = "character")
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("INFO.HGNC", "INFO.CQ", "POS", "CHROM", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "rauch"
    
    return(de_novos)
}

#' get de novo data for De Ligt et al. intellectual disability exome study
#' 
#' De novo mutation data sourced from supplementary table 3 from
#' De Ligt et al. (2012) N Engl J Med (2012) 367:1921-1929
#' doi: 10.1056/NEJMoa1206524
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
open_deligt_de_novos <- function() {
    
    de_novos = read.delim(file.path(DE_NOVO_DIR, "deligt_v2.txt"), header=TRUE, colClasses = "character")
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("INFO.HGNC", "INFO.CQ", "POS", "CHROM", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "deligt"
    
    return(de_novos)
}

#' get de novo data for the Epi4K epilepsy exome study
#' 
#' De novo mutation data sourced from supplementary table 2 from
#' Allen et al. (2013) Nature 501:217–221 
#' doi: 10.1038/nature12439
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
open_epi4k_de_novos <- function() {
    
    de_novos = read.delim(file.path(DE_NOVO_DIR, "epi4k_v2.txt"), header=TRUE, colClasses = "character")
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(de_novos$Type == "snv")
    de_novos$TYPE = "INDEL"
    de_novos$TYPE[TYPE.index] = "SNV"
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("INFO.HGNC", "INFO.CQ", "pos", "chrom", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "epi4k"
    
    return(de_novos)
}

#' get de novo data from autism exome studies
#' 
#' I think this data was obtained from studies published using the Simon's
#' Simplex Collection data (http://sfari.org/resources/simons-simplex-collection)
#' 
#' Supplementary table 2 (where the excel sheets for the probands and 
#' siblings have been combined) from:
#' Sanders et al. (2012) Nature 485:237-241 
#' doi: 10.1038/nature10945
#' 
#' Supplementary table 3 from:
#' O'Roak et al. (2012) Nature 485:246-250
#' doi: 10.1038/nature10989
#' 
#' Supplementary table 1 (where the non-coding SNVs have been excluded) from:
#' Iossifov et al. (2012) Neuron 74:285-299
#' doi: 10.1016/j.neuron.2012.04.009
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
open_autism_de_novos <- function() {
    
    de_novos = read.delim(file.path(DE_NOVO_DIR, "autism_v3_PJ.txt"), header=TRUE, colClasses = "character")
    
    # select only de novos in probands
    de_novos = de_novos[which(de_novos$pheno == "Pro"), ]
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(abs(nchar(de_novos$ref.1) - nchar(de_novos$var)) == 0)
    de_novos$TYPE = "INDEL"
    de_novos$TYPE[TYPE.index] = "SNV"
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("INFO.HGNC", "INFO.CQ", "pos", "CHROM", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "autism"
    
    return(de_novos)
}

#' get de novo data from Fromer et al. schizophrenia exome study
#' 
#' De novo mutation data sourced from Supplementary table 1:
#' Fromer et al. (2014) Nature 506:179–184
#' doi: 10.1038/nature12929
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
open_fromer_de_novos <- function() {
    
    de_novos = read.delim(file.path(DE_NOVO_DIR, "fromer_v2.txt"), header=TRUE, colClasses = "character")
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(abs(nchar(de_novos$Reference.allele) - nchar(de_novos$Alternate.allele)) == 0)
    de_novos$TYPE = "INDEL"
    de_novos$TYPE[TYPE.index] = "SNV"
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("INFO.HGNC", "INFO.CQ", "pos", "chrom", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "fromer"
    
    return(de_novos)
}

#' get de novo data from Zaidi et al. congenital heart disease exome study
#' 
#' De novo mutation data sourced from Supplementary table 4:
#' Zaidi et al. (2013) Nature 498:220–223
#' doi: 10.1038/nature12141
#' 
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
open_zaidi_de_novos <- function() {
    
    # could only include syndromic DNMs
    de_novos = read.delim(file.path(DE_NOVO_DIR, "zaidi_VEP.txt"), header=TRUE, colClasses = "character")
    
    # remove DNMs in controls
    de_novos = de_novos[-which(de_novos$Primary_Cardiac_Class == "Control"), ]
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(abs(nchar(de_novos$ref) - nchar(de_novos$alt)) != 0)
    TYPE.index = sort(unique(c(TYPE.index, which(de_novos$ref == "-"), which(de_novos$alt == "-"))))
    de_novos$TYPE = "SNV"
    de_novos$TYPE[TYPE.index] = "INDEL"
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("INFO.HGNC", "INFO.CQ", "pos", "chrom", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "zaidi"
    
    return(de_novos)
}
