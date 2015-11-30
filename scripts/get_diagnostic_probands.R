
library(argparse)
library(mupit)

get_options <- function() {
    parser = ArgumentParser()
    parser$add_argument("--ddd-1k-diagnoses", help="Path to DDD 1K diagnoses.",
        default="/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx")
    parser$add_argument("--de-novos", help="Path to DDD de novo dataset.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-10-12.txt")
    parser$add_argument("--low-pp-dnm", help="Path to low PP_DNM validations.",
        default="/nfs/users/nfs_j/jm33/de_novos.ddd_4k.validation_results.low_pp_dnm.2015-10-02.xlsx")
    parser$add_argument("--recessive-diagnoses", help="Path to additional recessive diagnoses.")
    parser$add_argument("--families", help="Path to families PED file.",
        default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt")
    parser$add_argument("--ddg2p", help="Path to DDG2P file.",
        default="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter")
    parser$add_argument("--out", help="Path to output file.",
        default=paste("/lustre/scratch113/projects/ddd/users/jm33/ddd_4k.diagnosed", Sys.Date(), "txt", sep="."))
        
    args = parser$parse_args()
    
    return(args)
}

get_ddd_diagnostic_cnvs <- function(path) {
    # load the CNVs
    reviewed_cnvs = gdata::read.xls(path, sheet="CNVs reviewed", stringsAsFactors=FALSE)
    
    # redefine the columns we want from the CNVs table
    cnvs = reviewed_cnvs[reviewed_cnvs$Diagnostic. > 0, ]
    cnvs$person_id = cnvs$DDDP_ID
    cnvs$chrom = cnvs$Chr
    cnvs$start_pos = as.numeric(cnvs$Start)
    cnvs$end_pos = as.numeric(cnvs$Stop)
    cnvs$ref_allele = NA
    cnvs$alt_allele = NA
    cnvs$hgnc = NA
    cnvs$inheritance = cnvs$Inheritance
    cnvs$type = "cnv"
    
    cnvs = cnvs[, c("person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "type")]
    
    return(cnvs)
}

get_ddd_diagnostic_snvs <- function(path) {
    # now load the diagnostic SNVs
    reviewed_snvs = gdata::read.xls(path, sheet="Exome variants reviewed", stringsAsFactors=FALSE)
    snvs = reviewed_snvs[reviewed_snvs$DECISION == "Yes", ]
    
    snvs$person_id = snvs$proband
    snvs$start_pos = as.numeric(snvs$position)
    alleles = strsplit(snvs$ref.alt_alleles, "/")
    snvs$ref_allele = sapply(alleles, "[", 1)
    snvs$alt_allele = sapply(alleles, "[", 2)
    snvs$end_pos = snvs$start_pos + nchar(snvs$ref_allele) - 1
    
    snvs$hgnc = snvs$gene
    snvs$inheritance = snvs$inheritance..DECIPHER.compatible.
    
    # determine the type of variant
    snvs$type = "snv"
    snvs$type[nchar(snvs$ref_allele) > 1 | nchar(snvs$alt_allele) > 1] = "indel"
    
    snvs = snvs[, c("person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "type")]
    
    return(snvs)
}

get_other_diagnostic_variants <- function(path) {
    # now load the diagnostic SNVs
    other = gdata::read.xls(path, sheet="UPD&Mosaicism", header=FALSE, stringsAsFactors=FALSE)
    names(other) = c("person_id", "type")
    other$chrom = NA
    other$start_pos = NA
    other$end_pos = NA
    other$ref_allele = NA
    other$alt_allele = NA
    other$hgnc = NA
    other$inheritance = NA
    
    other$type = gsub("Mosaicism", "mosaic_cnv", other$type)
    other$type = gsub("UPD", "uniparental_disomy", other$type)
    
    other = other[, c("person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "type")]
    
    return(other)
}

#' find diagnosed probands in the DDD study, to exclude them from our data
#'
#' @param path path to file defining diagnosed probands
#'
#' @export
#' @return A list containing vectors with DDD IDs, and sex of the diagnosed
#'     probands
get_previous <- function(path, families_path) {
    
    cnvs = get_ddd_diagnostic_cnvs(path)
    snvs = get_ddd_diagnostic_snvs(path)
    other = get_other_diagnostic_variants(path)
    
    diagnosed = rbind(cnvs, snvs, other)
    
    # relabel the inheritance types
    diagnosed$inheritance = gsub("DNM", "de_novo", diagnosed$inheritance)
    diagnosed$inheritance = gsub("De novo constitutive", "de_novo", diagnosed$inheritance)
    diagnosed$inheritance = gsub("Maternal", "maternal", diagnosed$inheritance)
    diagnosed$inheritance = gsub("Biparental", "biparental", diagnosed$inheritance)
    diagnosed$inheritance = gsub("[Pp]aternally inherited, constitutive in father", "paternal", diagnosed$inheritance)
    diagnosed$inheritance = gsub("[Mm]aternally inherited, constitutive in mother", "maternal", diagnosed$inheritance)
    diagnosed$inheritance = gsub("Biparental", "biparental", diagnosed$inheritance)
    
    # TODO: get the sex info for each proband
    families = read.table(families_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    families = families[, c("individual_id", "sex")]
    names(families) = c("person_id", "sex")
    
    families$sex = gsub("M", "male", families$sex)
    families$sex = gsub("F", "female", families$sex)
    
    diagnosed = merge(diagnosed, families, by="person_id", all.x=TRUE)
    
    return(diagnosed)
}

get_low_pp_dnm_validations <- function(path) {
    low_pp_dnm = gdata::read.xls(path, sheet="Summary_Final_forDB", stringsAsFactors=FALSE)
    
    # only select the most useful columns
    low_pp_dnm = low_pp_dnm[, c("ID", "CHR", "POS", "REF", "ALT", "TYPE",
        "manual_score")]
        
    names(low_pp_dnm) = c("person_id", "chrom", "start_pos", "ref_allele",
        "alt_allele", "type", "status")
    
    # recode the validation status codes to more understandable codes
    low_pp_dnm$status[low_pp_dnm$status == "dnm"] = "de_novo"
    low_pp_dnm$status[low_pp_dnm$status == "dnm_low_alt"] = "de_novo"
    low_pp_dnm$status[low_pp_dnm$status == "fp"] = "false_positive"
    low_pp_dnm$status[low_pp_dnm$status == "inherited_pat"] = "inherited"
    low_pp_dnm$status[low_pp_dnm$status == "p/u"] = "uncertain"
    low_pp_dnm$status[low_pp_dnm$status == "unclear"] = "uncertain"
    low_pp_dnm$status[low_pp_dnm$status == "parental_mosaic"] = "de_novo"
    
    low_pp_dnm$type[low_pp_dnm$type == "DENOVO-SNP"] = "snv"
    low_pp_dnm$type[low_pp_dnm$type == "DENOVO-INDEL"] = "indel"
    
    return(low_pp_dnm)
}

get_current_de_novos <- function(path) {
    # load the current set of de novos to be analysed
    variants = read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    # standardise the SNV or INDEL flag
    variants$type = "indel"
    variants$type[variants$var_type == "DENOVO-SNP"] = "snv"
    
    # standardise the columns, and column names
    variants$person_id = variants$person_stable_id
    variants$start_pos = variants$pos
    variants$ref_allele = variants$ref
    variants$alt_allele = variants$alt
    variants$end_pos = variants$start_pos + nchar(variants$ref_allele) - 1
    
    variants$hgnc = variants$symbol
    variants$inheritance = "de_novo"
    
    variants$sex = gsub("M", "male", variants$sex)
    variants$sex = gsub("F", "female", variants$sex)
    
    return(variants)
}

#' checks if a sites has a match in a previous dataset
#'
#' Some de novos in the current dataset were also present in a previous
#' datafreeze. We need to remove these so as to not double count diagnostic
#' variants. Unfortunately, som sites have shifted location slightly such as
#' indels, which can be difficult to locate.
#'
#'
check_for_match <- function(site, initial, pos=FALSE) {
    
    # set the missing match type dependent on whether we are looking for whether
    # we have a match, or the matching position.
    if (pos) {no_match = NA} else {no_match = FALSE}
    
    vars = initial[initial[["person_id"]] == site[["person_id"]], ]
    vars = vars[vars[["chrom"]] == site[["chrom"]] & !is.na(vars[["chrom"]]), ]
    
    if (nrow(vars) == 0) { return(no_match) }
    
    delta = abs(vars[["start_pos"]] - as.integer(site[["start_pos"]]))
    
    if (min(delta) > 20) { return(no_match) }
    
    close = delta == min(delta)
    
    if (sum(close) == 1) {
        if (!pos) { return(TRUE) } else { return(vars[["start_pos"]][close]) }
    } else { return(no_match) }
}

#' find probands likely to have diagnoses, to exclude them from our data
#'
#' We assume DDD samples that have
#'
#' @param path path to file defining diagnosed probands
#' @export
#'
#' @return A list containing vectors with DDD IDs, and sex of the diagnosed
#'     probands
get_ddd_diagnosed <- function(diagnosed_path, de_novo_path,
    low_pp_dnm_validations_path, ddg2p_path, families_path, recessive_path) {
    
    initial_diagnosed = get_previous(diagnosed_path, families_path)
    
    ddg2p = load_ddg2p(ddg2p_path)
    variants = get_current_de_novos(de_novo_path)
    
    # the candidates with low pp_dnm (< 0.9) in DDG2P genes were attempted to
    # validate. Those that did, we can swap the pp_dnm to 1, since these are now
    # high confidence de novos
    low_pp_dnm = get_low_pp_dnm_validations(low_pp_dnm_validations_path)
    variants = merge(variants, low_pp_dnm[, c("person_id", "chrom", "start_pos", "status")],
        by=c("person_id", "chrom", "start_pos"), all.x=TRUE)
    variants$pp_dnm[variants$status == "de_novo"] = 1
    
    # get the set of de novos from the current dataset that are likely to be
    # diagnostic. These are de novos in genes with dominant modes of inheritance,
    # or chrX de novos in males in genes with hemizygous mode of inheritance,
    # and where the site is high confidence (as determined by having a high
    # pp_dnm, or a missing pp_dnm)
    likely_diagnostic = variants[(variants$hgnc %in% ddg2p$gene[ddg2p$dominant] |
        (variants$hgnc %in% ddg2p$gene[ddg2p$hemizygous] & variants$sex == "male"))
        & (variants$pp_dnm > 0.1 | is.na(variants$pp_dnm)), ]
    likely_diagnostic = likely_diagnostic[, c("person_id", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "inheritance", "type",
        "sex")]
    
    # remove the sites from the likely diagnoses that form part of a confirmed
    # diagnosis
    in_prev = apply(likely_diagnostic, 1, check_for_match, initial=initial_diagnosed)
    likely_diagnostic = likely_diagnostic[!in_prev, ]
    
    diagnosed = rbind(initial_diagnosed, likely_diagnostic)
    
    if (!is.null(recessive_path)) {
        recessive = read.table(recessive_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        diagnosed = rbind(diagnosed, recessive)
    }
    
    return(diagnosed)
}

main <- function() {
    
    args = get_options()
    diagnosed = get_ddd_diagnosed(args$ddd_1k_diagnoses, args$de_novos,
        args$low_pp_dnm, args$ddg2p, args$families, args$recessive_diagnoses)
    
    write.table(diagnosed, file=args$out, sep="\t", quote=FALSE, row.names=FALSE)
}

main()
