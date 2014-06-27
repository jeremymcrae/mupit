# functions to open de novo data from different datasets, and standardise the
# format, so that datasets can be combined
# 
# 

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")
SRC_DIR = file.path(CODE_DIR, "src")
DE_NOVO_DIR = file.path(DATA_DIR, "de_novo_datasets")


open_ddd_denovos <- function(diagnosed, sample_ids=NA) {
    # read in DNMs, genes, CQ and type from file
    de_novos = read.delim(file.path(DE_NOVO_DIR, "DNG_Variants_20Feb2014_NonRed_Clean_NoTwins_NoClusters.txt"), header=T, colClasses = "character")
    
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

open_rauch_de_novos <- function() {
    de_novos = read.delim(file.path(DE_NOVO_DIR, "rauch_v2.txt"), header=T, colClasses = "character")
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("INFO.HGNC", "INFO.CQ", "POS", "CHROM", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "rauch"
    
    return(de_novos)
}

open_deligt_de_novos <- function() {
    de_novos = read.delim(file.path(DE_NOVO_DIR, "deligt_v2.txt"), header=T, colClasses = "character")
    
    # standardise the columns, and column names
    de_novos = subset(de_novos, select = c("INFO.HGNC", "INFO.CQ", "POS", "CHROM", "TYPE"))
    names(de_novos) = c("HGNC", "CQ", "POS", "CHROM", "TYPE")
    de_novos$STUDY = "deligt"
    
    return(de_novos)
}

open_epi4k_de_novos <- function() {
    ########### Epi4K dataset ############
    de_novos = read.delim(file.path(DE_NOVO_DIR, "epi4k_v2.txt"), header=T, colClasses = "character")
    
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

open_autism_de_novos <- function() {
    ########### Autism dataset ###############
    de_novos = read.delim(file.path(DE_NOVO_DIR, "autism_v3_PJ.txt"), header=T, colClasses = "character")
    
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

open_fromer_de_novos <- function() {
    ########### Schizophrenia dataset ###############
    de_novos = read.delim(file.path(DE_NOVO_DIR, "fromer_v2.txt"), header=T, colClasses = "character")
    
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

open_zaidi_de_novos <- function() {
    ########### CHD dataset ###############
    # could only include syndromic DNMs
    de_novos = read.delim(file.path(DE_NOVO_DIR, "zaidi_VEP.txt"), header=T, colClasses = "character")
    
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
