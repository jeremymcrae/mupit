# functions to standardise allele codes into ref/alt format as nonambiguous
# IUPAC bases

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

#' fix duplication allele codes
#'
#' fix the allele column for variants with allele columns that are structured
#' like "dup"
#'
#' @param variants data frame of variants
#' @param allele_column name of column containing allele information
#'
#' @export
#' @return a data frame with chrom, start_pos, end_pos and allele columns.
fix_dup_alleles <-function(variants, allele_column) {
    
    temp = grepl("dup", variants[[allele_column]])
    
    # if we don't have any of these variants, simply return the vartiant table,
    # as it will cause an error if this continues to run on zero variants
    if (!(any(temp))) { return(variants) }
    
    variants$start_pos[temp] = as.numeric(variants$start_pos[temp]) - 1
    variants$end_pos[temp] = as.numeric(variants$end_pos[temp]) - 1
    
    ref_alleles = apply(variants[temp, ], 1, get_sequence_in_region)
    variants[[allele_column]][temp] = paste(ref_alleles, "/", ref_alleles, ref_alleles, sep = "")
    
    variants$start_pos[temp] = as.numeric(variants$start_pos[temp]) + 1
    
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
