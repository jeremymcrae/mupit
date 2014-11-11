# functions to clean de novo data, and identify variant consequences

#' find the most severe VEP consequence for a variant
#' 
#' @param variant data frame or list for a variant, containing columns named
#'     "chrom", "start_pos", "end_pos", and "alt_allele" code for a single variant
#' @param build genome build to find consequences on
#' @param verbose flag indicating whether to print variants as they are checked
#' 
#' @export
#' @return a character string containing the most severe consequence, as per VEP
#'     annotation formats.
#' 
#' @examples
#' get_most_severe_vep_consequence(data.frame(chrom=c("1"), 
#'     start_pos=c("1000000"), end_pos=c("1000000"), alt_allele=c("A")))
#' get_most_severe_vep_consequence(list(chrom="1", start_pos="1000000", 
#'     end_pos="1000000", alt_allele="A"))
get_most_severe_vep_consequence <- function(variant, build="grch37", verbose=FALSE) {
    
    # only tolerate the grch37 and grch38 genome builds, since they are the only
    # genome builds supported by the Ensembl REST API
    allowed_builds = c("grch37", "grch38")
    stopifnot( tolower(build) %in% allowed_builds )
    
    if (verbose) {
        print(paste("chr", variant$chrom, ":", variant$start_pos, " ", 
            variant$alt_allele, sep=""))
    }
    
    base_url = "rest.ensembl.org/vep/human/region/"
    if (build == "grch37") {
        base_url = paste("grch37", base_url, sep = ".")
    }
    
    url = paste("http://", base_url, variant$chrom, ":", variant$start_pos, ":",
        variant$end_pos, "/", variant$alt_allele, 
        "?content-type=application/json", sep = "")
    
    json = rjson::fromJSON(file=url)
    
    return(json[[1]]$most_severe_consequence)
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
prepare_coordinates <- function(variants, coordinate_column, allele_column) {
    
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
    # variants$temp = NULL
    
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
prepare_coordinates_with_allele <- function(variants, coordinate_column, allele_column) {
    
    # get the start pos and chromosome
    variants$temp = strsplit(toupper(variants[[coordinate_column]]), ":")
    variants$chrom = sapply(variants$temp, "[[", 1)
    variants$chrom = gsub("CHR", "", variants$chrom)
    variants$start_pos = sapply(variants$temp, "[[", 2)
    
    # extract the end pos, and allele code
    variants$temp = strsplit(variants[[allele_column]], "/")
    variants$ref_allele = sapply(variants$temp, "[[", 1)
    variants$end_pos = as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1
    variants$end_pos = as.character(variants$end_pos)
    variants$alt_allele = sapply(variants$temp, "[[", 2)
    
    # and drop the temp column
    # variants$temp = NULL
    
    return(variants)
}
