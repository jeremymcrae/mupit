
library(mupit)
library(gdata)

DIAGNOSED_PATH = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx"
DE_NOVO_PATH = "/nfs/users/nfs_j/jm33/apps/mupit/data-raw/de_novo_datasets/de_novos.ddd_4k.ddd_only.2015-09-15.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt"
DDG2P_PATH = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"
OUTPUT_PATH = "/lustre/scratch113/projects/ddd/users/jm33/ddd_4k.diagnosed.2015-09-15.txt"

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
    other = gdata::read.xls(path, sheet="UPD&amp;Mosaicism", header=FALSE, stringsAsFactors=FALSE)
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

check_for_match <- function(site, initial) {
    
    vars = initial[initial[["person_id"]] == site[["person_id"]], ]
    vars = vars[vars[["chrom"]] == site[["chrom"]] & !is.na(vars[["chrom"]]), ]
    
    if (nrow(vars) == 0) { return(FALSE) }
    
    delta = abs(vars[["start_pos"]] - as.integer(site[["start_pos"]]))
    
    if (min(delta) > 20) { return(FALSE) }
    
    close = delta == min(delta)
    
    return(sum(close) == 1)
}

load_dominant_ddg2p <- function(path) {
    # load the current DDG2P dataset, which is missing a field from the header
    ddg2p = read.table(path, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
    
    # fix the header problem
    cols = ddg2p[1, ]
    cols[8] = "description"
    names(ddg2p) = as.character(cols)
    ddg2p = ddg2p[-1, ]
    
    # restrict outrselves to the high-confidence genes with a dominant mode of
    # inheritance
    ddg2p = ddg2p[ddg2p$type != "Possible DD Gene", ]
    dominant = ddg2p[ddg2p$mode %in% c("Monoallelic", "X-linked dominant"), ]
    
    return(dominant)
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
get_ddd_diagnosed <- function(diagnosed_path, de_novo_path, ddg2p_path, families_path) {
    
    initial_diagnosed = get_previous(diagnosed_path, families_path)
    
    dominant = load_dominant_ddg2p(ddg2p_path)
    variants = get_current_de_novos(de_novo_path)
    
    # get the set of de novos from the current dataset that are likely to be
    # diagnostic. These are de novos in genes with dominant modes of inheritance,
    # and where the site is high confidence (as determined by having a high
    # pp_dnm, or a missing pp_dnm)
    likely_diagnostic = variants[variants$hgnc %in% dominant$gene & (variants$pp_dnm > 0.1 | is.na(variants$pp_dnm)), ]
    likely_diagnostic = likely_diagnostic[, c("person_id", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "inheritance", "type",
        "sex")]
    
    # remove the sites from the likely diagnoses that form part of a confirmed
    # diagnosis
    in_prev = apply(likely_diagnostic, 1, check_for_match, initial=initial_diagnosed)
    likely_diagnostic = likely_diagnostic[!in_prev, ]
    
    diagnosed = rbind(initial_diagnosed, likely_diagnostic)
    
    return(diagnosed)
}

main <- function() {
    diagnosed = get_ddd_diagnosed(DIAGNOSED_PATH, DE_NOVO_PATH, DDG2P_PATH, FAMILIES_PATH)
    
    write.table(diagnosed, file=OUTPUT_PATH, sep="\t", quote=FALSE, row.names=FALSE)
}

main()
