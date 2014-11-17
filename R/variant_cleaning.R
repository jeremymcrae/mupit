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
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], " ", 
            variant[["alt_allele"]], "    ", url, sep=""))
    }
    
    json = try(rjson::fromJSON(file=url), silent=TRUE)
    if (class(json) == "try-error") {json = rjson::fromJSON(file=ref_url)}
    
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
    
    # make sure the end position is suitable for the Ensembl REST API request
    if (as.numeric(variant[["end_pos"]]) < as.numeric(variant[["start_pos"]])) {
        variant[["end_pos"]] = variant[["start_pos"]]
    }
    
    # define parts of the URL
    url_start = paste("http://", base_url, variant[["chrom"]], ":", 
        variant[["start_pos"]], ":", variant[["end_pos"]], "/", sep="")
    content = "?feature=gene;content-type=application/json"
    url = paste(url_start, content, sep = "")
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], 
            "    ", url, sep=""))
    }
    
    json = rjson::fromJSON(file=url)
    
    # return blank string for variants not in genes
    if (length(json) > 0) {
        return(json[[1]]$external_name)
    }
    
    return("")
}

#' find genomic sequence within a region
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
get_sequence_in_region <- function(variant, build="grch37", verbose=FALSE) {
    
    # only tolerate the grch37 and grch38 genome builds, since they are the only
    # genome builds supported by the Ensembl REST API
    allowed_builds = c("grch37", "grch38")
    stopifnot( tolower(build) %in% allowed_builds )
    
    base_url = "rest.ensembl.org/sequence/region/human/"
    if (build == "grch37") {
        base_url = paste("grch37", base_url, sep = ".")
    }
    
    # make sure the end position is suitable for the Ensembl REST API request
    if (as.numeric(variant[["end_pos"]]) < as.numeric(variant[["start_pos"]])) {
        variant[["end_pos"]] = variant[["start_pos"]]
    }
    
    # define the URL
    url = paste("http://", base_url, variant[["chrom"]], ":", 
        variant[["start_pos"]], ":", variant[["end_pos"]], 
        ":1?content-type=application/json", sep="")
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], 
            "    ", url, sep=""))
    }
    
    json = rjson::fromJSON(file=url)
    
    # return blank string for variants not in genes
    if (length(json) > 0) {
        return(json$seq)
    }
    
    return("")
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

#' fix substitution allele codes
#' 
#' fix the allele column for variants with allele columns that are structured
#' like "sub(G-&gt;T)" (the "-&gt;" translates to ">", and sub(G->T), this is 
#' an excel to R conversion issue. Perhaps it is unicode from excel?
#' 
#' @param variants data frame of variants
#' @param allele_column name of column containing allele information
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
fix_sub_alleles <-function(variants, allele_column) {
    
    temp = grepl("sub", variants[[allele_column]])
    variants[[allele_column]][temp] = gsub("-&gt;|>", "/", variants[[allele_column]][temp])
    variants[[allele_column]][temp] = gsub("\\(|\\)", "", variants[[allele_column]][temp])
    variants[[allele_column]][temp] = gsub("sub", "", variants[[allele_column]][temp])
    
    return(variants)
}

#' fix deletion allele codes
#' 
#' fix the allele column for variants with allele columns that are structured
#' like "del(1)"
#' 
#' @param variants data frame of variants
#' @param allele_column name of column containing allele information
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
fix_del_alleles <- function(variants, allele_column) {
    
    temp = grepl("del", variants[[allele_column]])
    
    # if we don't have any of these variants, simply return the vartiant table,
    # as it will cause an error if this continues to run on zero variants
    if (!(any(temp))) { return(variants) }
    
    variants[[allele_column]][temp] = gsub("\\(|\\)", "", variants[[allele_column]][temp])
    variants[[allele_column]][temp] = gsub("del", "", variants[[allele_column]][temp])
    temp_distance = variants[[allele_column]][temp]
    
    # find the reference sequence at the site
    variants[[allele_column]][temp] = apply(variants[temp, ], 1, get_sequence_in_region)
    
    # find the sequence at the site + the distance of the deletion
    variants$end_pos[temp] = as.numeric(variants$end_pos[temp]) + as.numeric(temp_distance)
    variants[[allele_column]][temp] = paste(apply(variants[temp, ], 1, get_sequence_in_region), "/", variants[[allele_column]][temp],  sep="")
    
    return(variants)
}

#' fix insertion allele codes
#' 
#' fix the allele column for variants with allele columns that are structured
#' like "ins(1)"
#' 
#' @param variants data frame of variants
#' @param allele_column name of column containing allele information
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
fix_ins_alleles <-function(variants, allele_column) {
    
    temp = grepl("ins", variants[[allele_column]])
    
    # if we don't have any of these variants, simply return the vartiant table,
    # as it will cause an error if this continues to run on zero variants
    if (!(any(temp))) { return(variants) }
    
    variants[[allele_column]][temp] = gsub("\\(|\\)", "", variants[[allele_column]][temp])
    variants[[allele_column]][temp] = gsub("ins", "", variants[[allele_column]][temp])
    
    variants$start_pos[temp] = as.numeric(variants$start_pos[temp]) - 1
    variants$end_pos[temp] = as.numeric(variants$end_pos[temp]) - 1
    
    ref_alleles = apply(variants[temp, ], 1, get_sequence_in_region)
    variants[[allele_column]][temp] = paste(ref_alleles, "/", ref_alleles, variants[[allele_column]][temp], sep = "")
    
    variants$start_pos[temp] = as.numeric(variants$start_pos[temp]) + 1
    
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
#' prepare_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("del(1)")), "coord", "allele")
#' prepare_coordinates_with_allele(data.frame(coord=c("chr1:10000"), 
#'     allele=c("ins(ATG)")), "coord", "allele")
prepare_coordinates_with_allele <- function(variants, coordinate_column, allele_column) {
    
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
#' prepare_coordinates_with_hgvs_genomic(
#'    data.frame(hgvs=c("chr7:g.155556643G>A")), "hgvs")
#' prepare_coordinates_with_hgvs_genomic(data.frame(
#'    hgvs=c("chr3:g.11060365_11060365del")), "hgvs")
#' prepare_coordinates_with_hgvs_genomic(data.frame(
#'    hgvs=c("chr13:g.50057690_50057691insA")), "hgvs")
prepare_coordinates_with_hgvs_genomic <- function(variants, hgvs_column) {
    
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
#' @param coordinate_column name of column containing chromosome information
#' @param allele_column name of column containing allele information
#' 
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
#' 
#' @examples
#' prepare_zaiidi_coordinates(data.frame(chrom=c("1"), start_pos=c("100000"),
#'     alleles=c("A/G")), "alleles")
#' prepare_zaiidi_coordinates(data.frame(chrom=c("1"), start_pos=c("100000"),
#'     alleles=c("-AAAA")), "alleles")
#' prepare_zaiidi_coordinates(data.frame(chrom=c("1"), start_pos=c("100000"),
#'     alleles=c("+AAAA")), "alleles")
prepare_zaiidi_coordinates <- function(variants, allele_column) {
    
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
