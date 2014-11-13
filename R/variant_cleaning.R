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
#'     start_pos=c("1000000"), end_pos=c("1000000"), alt_allele=c("A"), 
#'     ref_allele=c("G")))
#' get_most_severe_vep_consequence(list(chrom="1", start_pos="1000000", 
#'     end_pos="1000000", alt_allele="A", ref_allele="G"))
get_most_severe_vep_consequence <- function(variant, build="grch37", verbose=FALSE) {
    
    # only tolerate the grch37 and grch38 genome builds, since they are the only
    # genome builds supported by the Ensembl REST API
    allowed_builds = c("grch37", "grch38")
    stopifnot( tolower(build) %in% allowed_builds )
    
    # get the correct allele for deletions for the Ensembl REST API requests
    if (variant[["alt_allele"]] == "") {variant[["alt_allele"]] = "-"}
    
    base_url = "rest.ensembl.org/vep/human/region/"
    if (build == "grch37") {
        base_url = paste("grch37", base_url, sep = ".")
    }
    
    # define parts of the URL
    url_start = paste("http://", base_url, variant[["chrom"]], ":", 
        variant[["start_pos"]], ":", variant[["end_pos"]], "/", sep="")
    content = "?content-type=application/json"
    
    # define the URLs for the ref and alt alleles
    url = paste(url_start, variant[["alt_allele"]], content, sep = "")
    ref_url = paste(url_start, variant[["ref_allele"]], content, sep = "")
    
    json = try(rjson::fromJSON(file=url), silent=TRUE)
    if (class(json) == "try-error") {json = rjson::fromJSON(file=ref_url)}
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], " ", 
            variant[["alt_allele"]], "    ", url, sep=""))
    }
    
    return(json[[1]]$most_severe_consequence)
}

#' find the hgnc symbol overlapping a variant position
#' 
#' @param variant data frame or list for a variant, containing columns named
#'     "chrom", "start_pos", and "end_pos" for a single variant
#' @param build genome build to find consequences on
#' @param verbose flag indicating whether to print variants as they are checked
#' 
#' @export
#' @return a character string containing the HGNC symbol.
#' 
#' @examples
#' get_gene_id_for_variant(data.frame(chrom=c("1"), start_pos=c("1000000"),
#'     end_pos=c("1000000")))
#' get_gene_id_for_variant(list(chrom="1", start_pos="1000000", 
#'     end_pos="1000000"))
get_gene_id_for_variant <- function(variant, build="grch37", verbose=FALSE) {
    
    # only tolerate the grch37 and grch38 genome builds, since they are the only
    # genome builds supported by the Ensembl REST API
    allowed_builds = c("grch37", "grch38")
    stopifnot( tolower(build) %in% allowed_builds )
    
    base_url = "rest.ensembl.org/overlap/region/human/"
    if (build == "grch37") {
        base_url = paste("grch37", base_url, sep = ".")
    }
    
    # define parts of the URL
    url_start = paste("http://", base_url, variant[["chrom"]], ":", 
        variant[["start_pos"]], ":", variant[["end_pos"]], "/", sep="")
    content = "?feature=gene;content-type=application/json"
    url = paste(url_start, content, sep = "")
    
    json = rjson::fromJSON(file=url)
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], 
            "    ", url, sep=""))
    }
    
    return(json[[1]]$external_name)
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
#' prepare_coordinates(data.frame(coord=c("chr1:10000:A"), 
#'     allele=c("chr1:10000:A:G")), "coord", "allele")
#' prepare_coordinates(data.frame(coord=c("chr1:10000:A"), 
#'     allele=c("chr1:10003:ATGC:G")), "coord", "allele")
prepare_coordinates <- function(variants, coordinate_column, allele_column) {
    
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
#' prepare_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("A/G")), "coord", "allele")
#' prepare_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("ATGC/G")), "coord", "allele")
#' prepare_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("sub(G-&gt;T)")), "coord", "allele")
prepare_coordinates_with_allele <- function(variants, coordinate_column, allele_column) {
    
    # make sure the required columns are as character
    variants[[coordinate_column]] = as.character(variants[[coordinate_column]])
    variants[[allele_column]] = as.character(variants[[allele_column]])
    
    # get the start pos and chromosome
    variants$temp = strsplit(toupper(variants[[coordinate_column]]), ":")
    variants$chrom = sapply(variants$temp, "[[", 1)
    variants$chrom = gsub("CHR", "", variants$chrom)
    variants$start_pos = sapply(variants$temp, "[[", 2)
    
    # fix the allele column for variants with allele columns that are structured
    # like "sub(G-&gt;T)" (the "-&gt;" translates to ">", and sub(G->T), this is 
    # an excel to R conversion issue. Perhaps it is unicode from excel?
    temp = grepl("sub", variants[[allele_column]])
    variants[[allele_column]][temp] = gsub("-&gt;", "/", variants[[allele_column]][temp])
    variants[[allele_column]][temp] = gsub("\\(|\\)", "", variants[[allele_column]][temp])
    variants[[allele_column]][temp] = gsub("sub", "", variants[[allele_column]][temp])
    
    # extract the end pos, and allele code
    variants$temp = strsplit(variants[[allele_column]], "/")
    variants$ref_allele = sapply(variants$temp, "[[", 1)
    variants$end_pos = as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1
    variants$end_pos = as.character(variants$end_pos)
    variants$alt_allele = sapply(variants$temp, "[[", 2)
    
    # and drop the temp column
    variants$temp = NULL
    
    return(variants)
}

#' correct alt alleles which have been encoded as an IUPAC ambiguous base
#' 
#' Sometimes the de novo variants from studies provide alt alleles as "R", or
#' "Y", which indicate ambigous bases. We can identify the correct alt base by
#' comparison with the reference allele.
#' 
#' @param variants data frame of variants
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
#' 
#' @examples
#' prepare_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("A/G")), "coord", "allele")
#' prepare_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("ATGC/G")), "coord", "allele")
fix_het_alleles <- function(variants) {
    iupac = list("R" = c("A", "G"), "Y" = c("C", "T"), "S" = c("G", "C"),
                 "W" = c("A", "T"), "K" = c("G", "T"), "M" = c("A", "C"))
    
    # we need to correct
    for (row_num in 1:nrow(variants)) {
        alt = variants$alt_allele[row_num]
        ref = variants$ref_allele[row_num]
        
        # figure the correct base from the ambigous base which is not the 
        # reference allele.
        if (alt %in% names(iupac)) {
            alt = setdiff(iupac[[alt]], ref)
            variants$alt_allele[row_num] = alt
        }
    }
    
    return(variants)
}
