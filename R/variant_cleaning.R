# functions to clean de novo data, and identify variant consequences

#' extract genomic coordinates for variants from a dataframe
#' 
#' Genomic coordinates can be encoded in a variety of ways, this function
#' cleans up coordinates where one column
#' 
#' @param variants data frame of variants
#' @param coordinate_column name of column containing chromosome information
#' @param allele_column name of column containing allele information
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns. 
#' 
#' @examples
#' fix_coordinates(data.frame(coord=c("chr1:10000:A"), 
#'     allele=c("chr1:10000:A:G")), "coord", "allele")
#' fix_coordinates(data.frame(coord=c("chr1:10000:A"), 
#'     allele=c("chr1:10003:ATGC:G")), "coord", "allele")
fix_coordinates <- function(variants, coordinate_column, allele_column) {
    
    # make sure the required columns are as character
    variants[[coordinate_column]] = as.character(variants[[coordinate_column]])
    variants[[allele_column]] = as.character(variants[[allele_column]])
    
    # get the start pos and chromosome
    variants$temp = strsplit(toupper(variants[[coordinate_column]]), ":")
    variants$chrom = sapply(variants$temp, "[[", 1)
    variants$chrom = gsub("CHR", "", variants$chrom)
    variants$start_pos = sapply(variants$temp, "[[", 2)
    
    # extract the end pos, and allele code
    variants$temp = strsplit(variants[[allele_column]], ":")
    variants$end_pos = sapply(variants$temp, "[[", 2)
    variants$ref_allele = sapply(variants$temp, "[[", 3)
    variants$alt_allele = sapply(variants$temp, "[[", 4)
    
    # and drop the temp column
    variants$temp = NULL
    
    return(variants)
}

#' extract genomic coordinates for variants from a dataframe
#' 
#' Genomic coordinates can be encoded in a variety of ways, this function
#' cleans up coordinates where one column
#' 
#' @param variants data frame of variants
#' @param coordinate_column name of column containing chromosome information
#' @param allele_column name of column containing allele information
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
#' 
#' @examples
#' fix_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("A/G")), "coord", "allele")
#' fix_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("ATGC/G")), "coord", "allele")
#' fix_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("sub(G-&gt;T)")), "coord", "allele")
#' fix_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("del(1)")), "coord", "allele")
#' fix_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("ins(ATG)")), "coord", "allele")
fix_coordinates_with_allele <- function(variants, coordinate_column, allele_column) {
    
    # make sure the required columns are as character
    variants[[coordinate_column]] = as.character(variants[[coordinate_column]])
    variants[[allele_column]] = as.character(variants[[allele_column]])
    
    # get the start pos and chromosome
    variants$temp = strsplit(toupper(variants[[coordinate_column]]), ":")
    variants$chrom = sapply(variants$temp, "[[", 1)
    variants$chrom = gsub("CHR", "", variants$chrom)
    variants$start_pos = sapply(variants$temp, "[[", 2)
    variants$end_pos = variants$start_pos
    
    variants = fix_sub_alleles(variants, allele_column)
    variants = fix_del_alleles(variants, allele_column)
    variants = fix_ins_alleles(variants, allele_column)
    
    # extract the end pos, and allele code
    variants$temp = strsplit(variants[[allele_column]], "/")
    variants$ref_allele = sapply(variants$temp, "[[", 1)
    variants$end_pos = as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1
    variants$end_pos = as.character(variants$end_pos)
    variants$alt_allele = sapply(variants$temp, "[[", 2)
    
    # drop the temp column
    variants$temp = NULL
    
    return(variants)
}

#' extract genomic coordinates for variants encoded as HGVS genomic
#' 
#' 
#' @param variants data frame of variants
#' @param hgvs_column name of column containing HGVS coordinate
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
#' 
#' @examples
#' fix_coordinates_with_hgvs_genomic(
#'    data.frame(hgvs=c("chr7:g.155556643G>A")), "hgvs")
#' fix_coordinates_with_hgvs_genomic(data.frame(
#'    hgvs=c("chr3:g.11060365_11060365del")), "hgvs")
#' fix_coordinates_with_hgvs_genomic(data.frame(
#'    hgvs=c("chr13:g.50057690_50057691insA")), "hgvs")
fix_coordinates_with_hgvs_genomic <- function(variants, hgvs_column) {
    
    # make sure the required columns are as character
    variants[[hgvs_column]] = as.character(variants[[hgvs_column]])
    
    # get the chromosome
    variants$temp = strsplit(variants[[hgvs_column]], ":")
    variants$chrom = sapply(variants$temp, "[[", 1)
    variants$chrom = gsub("chr", "", variants$chrom)
    
    # prepare the position and allele strings
    variants$temp = sapply(variants$temp, "[[", 2)
    variants$temp = strsplit(variants$temp, "\\.")
    variants$temp = sapply(variants$temp, "[[", 2)
    
    # get the position strings
    variants$start_pos = strsplit(variants$temp, "[a-zA-Z]")
    variants$start_pos = sapply(variants$start_pos, "[[", 1)
    variants$end_pos = variants$start_pos
    indels = grepl("_", variants$start_pos)
    variants$start_pos[indels] = sapply(strsplit(variants$start_pos[indels], "_"), "[[", 1)
    variants$end_pos[indels] = sapply(strsplit(variants$end_pos[indels], "_"), "[[", 2)
    
    # get the allele string
    variants$temp_alleles = strsplit(variants$temp, "[0-9_]+")
    variants$temp_alleles = sapply(variants$temp_alleles, "[[", 2)
    
    # fix the substitution
    sub = grepl("del[a-zA-Z0-9]+ins", variants$temp_alleles)
    variants$temp_alleles[sub] = gsub("del", "", variants$temp_alleles[sub])
    variants$temp_alleles[sub] = gsub("ins", "/", variants$temp_alleles[sub])
    
    # deletions
    del = grepl("del", variants$temp_alleles)
    dist = as.numeric(variants$end_pos[del]) - as.numeric(variants$start_pos[del]) + 1
    variants$end_pos[del] = variants$start_pos[del]
    variants$temp_alleles[del] = paste("del", "(", dist, ")", sep = "")
    variants = fix_del_alleles(variants, "temp_alleles")
    
    # insertions
    ins = grepl("ins", variants$temp_alleles)
    variants$end_pos[ins] = variants$start_pos[ins]
    variants$temp_alleles[ins] = paste("ins", "(", gsub("ins", "", variants$temp_alleles[ins]), ")", sep = "")
    variants = fix_ins_alleles(variants, "temp_alleles")
    
    # duplications
    dup = grepl("dup", variants$temp_alleles)
    variants = fix_dup_alleles(variants, "temp_alleles")
    
    # standardise the allele separator
    variants$temp_alleles = gsub(">", "/", variants$temp_alleles)
    
    # get the ref and alt alleles
    variants$temp_alleles = strsplit(variants$temp_alleles, "/")
    variants$ref_allele = sapply(variants$temp_alleles, "[[", 1)
    variants$end_pos = as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1
    variants$end_pos = as.character(variants$end_pos)
    variants$alt_allele = sapply(variants$temp_alleles, "[[", 2)
    
    # drop the temp columns
    variants$temp = NULL
    variants$temp_alleles = NULL
    
    return(variants)
}

#' extract genomic coordinates for variants from a dataframe
#' 
#' Genomic coordinates can be encoded in a variety of ways, this function
#' cleans up coordinates where one column
#' 
#' @param variants data frame of variants
#' @param allele_column name of column containing allele information
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
#' 
#' @examples
#' fix_zaiidi_coordinates(data.frame(chrom=c("1"), start_pos=c("100000"),
#'     alleles=c("A/G")), "alleles")
#' fix_zaiidi_coordinates(data.frame(chrom=c("1"), start_pos=c("100000"),
#'     alleles=c("-AAAA")), "alleles")
#' fix_zaiidi_coordinates(data.frame(chrom=c("1"), start_pos=c("100000"),
#'     alleles=c("+AAAA")), "alleles")
fix_zaiidi_coordinates <- function(variants, allele_column) {
    
    # make sure the required columns are as character
    variants[[allele_column]] = as.character(variants[[allele_column]])
    variants$start_pos = as.character(variants$start_pos)
    variants$end_pos = variants$start_pos
    
    # fix the one indel
    indels = grepl("-[a-zA-Z0-9\\/]+\\+", variants[[allele_column]])
    variants[[allele_column]][indels] = gsub("-|\\+", "", variants[[allele_column]][indels])
    
    # fix the deletions
    dels = grepl("-", variants[[allele_column]])
    variants[[allele_column]][dels] = gsub("-", "", variants[[allele_column]][dels])
    variants[[allele_column]][dels] = paste("del(", nchar(variants[[allele_column]][dels]), ")", sep="")
    variants = fix_del_alleles(variants, allele_column)
    
    # fix the insertions
    ins = grepl("\\+", variants[[allele_column]])
    variants[[allele_column]][ins] = gsub("\\+", "", variants[[allele_column]][ins])
    variants[[allele_column]][ins] = paste("ins(", variants[[allele_column]][ins], ")", sep="")
    variants = fix_ins_alleles(variants, allele_column)
    
    # get the ref and alt alleles
    variants[[allele_column]] = gsub(">", "/", variants[[allele_column]])
    variants[[allele_column]] = strsplit(variants[[allele_column]], "/")
    variants$ref_allele = sapply(variants[[allele_column]], "[[", 1)
    variants$end_pos = as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1
    variants$end_pos = as.character(variants$end_pos)
    variants$alt_allele = sapply(variants[[allele_column]], "[[", 2)
    
    variants[[allele_column]] = NULL
    
    return(variants)
}
