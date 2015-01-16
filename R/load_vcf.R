# functions to help loading VCFs

#' open all the variants found within a single-sample VCF
#' 
#' @param path path to VCF
#' @export
#' 
#' @return dataframe of all variants in the VCF, modified so that the sample ID
#'     is an additional column
load_vcf <- function(path) {
    vcf = readLines(path)
    
    # strip the header lines, convert the final header to a normal header line
    vcf = vcf[!grepl("^##", vcf)]
    vcf = gsub("#", "", vcf)
    
    # load the variants as a dataframe
    vcf = read.table(textConnection(vcf), sep="\t", header=TRUE, 
        colClasses="character")
    
    if (ncol(vcf) != 10) { stop("only opens single-sample VCFs, check columns") }
    
    # get the sample id
    vcf$sample_id = character(nrow(vcf))
    
    if (nrow(vcf) > 0) {
        vcf$sample_id = names(vcf)[10]
    }
    
    # standardise the format values column, away from being the sample ID
    names(vcf)[10] = "format_values"
    
    return(vcf)
}

#' open all the variants found within a set of single-sample VCFs
#' 
#' @param vcf_paths vector of paths to VCFs
#' @export
#' 
#' @return dataframe of all variants in the VCFs
get_variants <- function(vcf_paths) {
    variants = vector("list", length=length(vcf_paths))
    
    # report progress when loading the VCFs, ~70s for 4300 probands.
    cat("loading VCFs\n")
    pb = txtProgressBar(min=0, max=length(vcf_paths), style=3)
    
    # open each VCF, and add the variants to a list entry
    for (pos in 1:length(vcf_paths)) {
        if (pos %% 20 == 0) {setTxtProgressBar(pb, pos)}
        
        path = vcf_paths[pos]
        variants[[pos]] = load_vcf(path)
    }
    
    close(pb)
    
    variants = dplyr::rbind_all(variants)
    
    return(variants)
}

#' get HGNC symbols from the INFO field of VCF variants
#' 
#' @param info vector of VCF INFO entries for different variants
#' @export
#' 
#' @return vector of HGNC symbols, in same order as the INFO vector
get_hgnc_from_vcf_info <- function(info) {
    
    if (length(info) == 0) { return(rep(NA, length(info))) }
        
    # sometimes the variants have HGNC info fields, while other times they
    # lack HGNC, but do have HGNC_ALL. We prioritise the HGNC field.
    hgnc_info = rep(NA, length(info))
    hgnc_info[grepl("HGNC_ALL=", info)] = "HGNC_ALL"
    hgnc_info[grepl("HGNC=", info)] = "HGNC"
    hgnc_pattern = paste(hgnc_info, "=[a-zA-Z0-9&-]+;", sep="")
    
    # get the start and end positions of the HGNC in the INFO field
    hgnc_pos = stringr::str_locate_all(info, hgnc_pattern)
    temp = data.frame("start"=numeric(0), end=numeric(0))
    hgnc_pos = data.frame(t(sapply(hgnc_pos, rbind, temp)))
    
    # trim down to the HGNC symbol
    start_pos = unlist(hgnc_pos$start) + nchar("HGNC=")
    hgnc_all = hgnc_info == "HGNC_ALL"
    start_pos[hgnc_all] = unlist(hgnc_pos$start[hgnc_all]) + nchar("HGNC_ALL=")
    end_pos = unlist(hgnc_pos$end) - 1
    
    hgnc = substr(info, start_pos, end_pos)
    
    return(hgnc)
}
