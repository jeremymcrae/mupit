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
    variants = as.data.frame(variants)
    
    return(variants)
}
