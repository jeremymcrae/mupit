### script to identify de novo mutations for validation

library(mupit)

# get the functions used for de novo filtering
options(run.main=FALSE)
source("/nfs/users/nfs_j/jm33/apps/denovo_filter/filter_DNG.R")

ALL_PROBANDS_RESULTS = "/nfs/users/nfs_j/jm33/apps/mupit/results/de_novos.ddd_4k.with_diagnosed.all.txt"
ALL_PROBANDS_RESULTS_WITHOUT_DIAGNOSED = "/nfs/users/nfs_j/jm33/apps/mupit/results/de_novos.ddd_4k.without_diagnosed.all.txt"
DECIPHER_ID_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/person_sanger_decipher.txt"

DE_NOVOS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/denovo_gear_trios_extracted_passed_variants_18.12.14.tsv"
VALIDATION_PATH = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/DNG_Validation_1133trios_20140130.tsv"


#' loads the results of significane testing from gene tests
#'
#' @param with_diagnosed_path path to file containing results from testing for
#'     associations, using the dataset that contains the probands with diagnoses
#' @param without_diagnosed_path path to file containing results from testing for
#'     associations, using the dataset that contains the probands without diagnoses
load_genes <- function(with_diagnosed_path, without_diagnosed_path) {
    
    with_diagnosed = read.table(with_diagnosed_path, header=TRUE)
    
    # indicate which genes achieve genome-wide significance
    with_diagnosed$significant = with_diagnosed$p_min * 38500 < 0.05
    
    # restrict it to the minimum columns
    with_diagnosed = with_diagnosed[, c("hgnc", "in_ddg2p", "meta.lof_indel",
        "meta.lof_snv", "meta.missense_indel", "meta.missense_snv", "p_min")]
    
    # now open the P values for the tests without the diagnosed probands
    without_diagnosed = read.table(without_diagnosed_path, header=TRUE)
    without_diagnosed = without_diagnosed[, c("hgnc", "p_min")]
    names(without_diagnosed) = c("hgnc", "p_min_without_diagnosed")
    
    genes = merge(with_diagnosed, without_diagnosed, by="hgnc")
    names(genes) = gsub("meta.", "", names(genes))
    genes = genes[order(genes$p_min), ]
    
    return(genes)
}

#' show the details of the de novos for a gene
#'
#' @param hgnc HGNC symbol for a gene
#' @param de novos dataframe of fully annotated denovogear output
#' @param validated dataframe of variants that have had validation efforts
examine_de_novos <- function(hgnc, de_novos, validated) {
    
    # select the rows corresponding to the gene
    vars = de_novos[de_novos$symbol == hgnc, ]
    
    validated = validated[, c("person_stable_id", "chr", "pos", "validation_result")]
    names(validated) = c("person_stable_id", "chrom", "pos", "validated")
    
    vars = merge(vars, validated, by=c("person_stable_id", "chrom", "pos"), all.x=TRUE)
    
    # and only include certain columns
    vars = vars[, c("symbol", "person_stable_id", "gender", "chrom", "pos",
    "ref", "alt", "consequence", "max_af", "pp_dnm", "rd_child", "rd_dad",
    "rd_mom", "child_alt_prp", "max_alt_in_parent", "child.REF.F",
    "child.REF.R", "child.ALT.F", "child.ALT.R", "father.REF.F",
    "father.REF.R", "father.ALT.F", "father.ALT.R", "mother.REF.F",
    "mother.REF.R", "mother.ALT.F", "mother.ALT.R", "REF.F", "REF.R",
    "ALT.F", "ALT.R", "parent.ALT", "parent.REF", "PA_pval_site", "SB_pval",
    "PA_pval_gene", "min.parent.ALT","gene.ALT", "gene.REF", "overall.pass",
    "validated")]
    
    print(vars)
}

#' count the numner of de novos in a given gene in the DDD dataset
#'
#' @param hgnc HGNC symbol for gene
#' @param functional whether to exclude non-functional (ie synonymous and UTR
#'     variants)
#' @param validated dataset of validated de novos, if used we exclude the
#'     validated de novos from the final count
#' @param return_data boolean for whether to return the full dataframe that
#'     includes the
count_in_ddd <- function(hgnc, functional=TRUE, validated=NULL, return_data=FALSE) {
    
    vars = mupit::ddd_de_novos[mupit::ddd_de_novos$hgnc %in% hgnc, ]
    
    if (functional) { vars = vars[vars$consequence != "synonymous_variant", ] }
    
    if (!is.null(validated)) {
        validated = validated[, c("person_stable_id", "chr", "pos", "validation_result")]
        names(validated) = c("person_id", "chrom", "start_pos", "validated")
        
        vars = merge(vars, validated, by=c("person_id", "chrom", "start_pos"), all.x=TRUE)
        vars = vars[is.na(vars$validated), ]
    }
    
    if (!return_data) {
        return(nrow(vars))
    } else {
        return(vars)
    }
}

#' load the full set of denovogear calls, and annotate with filtering status
#'
#' @param path path to the annotated denovogear output
#'
#' @return dataframe of de novo variants, including filtering status
get_raw_de_novos <- function(path) {
    de_novos = read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    de_novos = fix_maf(de_novos)
    de_novos = fix_missing_gene_symbols(de_novos)
    de_novos = extract_alt_and_ref_counts(de_novos)
    de_novos = count_site_and_gene_recurrence(de_novos)
    
    site_results = test_sites(de_novos)
    de_novos = merge(de_novos, site_results, by="key", all.x=TRUE)
    
    # get the minimum alternate allele count from the parents
    alts = data.frame(de_novos$mother.ALT.F + de_novos$mother.ALT.R, de_novos$father.ALT.F + de_novos$father.ALT.R)
    de_novos$min.parent.ALT = apply(alts, 1, min)
    
    # test whether any genes have an excess of parental alts
    gene_results = test_genes(de_novos)
    de_novos = merge(de_novos, gene_results, by="symbol", all.x=TRUE)
    
    de_novos = set_filter_flag(de_novos, keep_all=TRUE)
    
    return(de_novos)
}

#' find the variants in the potentially novel genes to validate
#'
#' @param genes dataframe of results from testing for de novo enrichment,
#'     clustering, and HPO similarity.
#' @param validations dataframe of validation results from the 1133 trios
get_variants_to_validate <- function(genes, validations, path) {
    
    # restrict the gene results to genes which are not known developmental
    # disorders, (ie are potentially novel genes)
    genes = genes[!genes$in_ddg2p, ]
    
    # restrict the validation dataframe to the minimum required to identify the
    # variant, and the validation result
    valid = validations[, c("person_stable_id", "chr", "pos", "validation_result")]
    names(valid) = c("person_id", "chrom", "start_pos", "validated")
    
    temp = merge(mupit::ddd_de_novos, valid, by=c("person_id", "chrom", "start_pos"), all.x=TRUE)
    
    # find the top 100 genes (a reasonable set to perform validations on),
    # then remove variants which have previously been checked (either giving a
    # validated or invalidated result)
    temp = temp[temp$hgnc %in% genes$hgnc[1:100], ]
    to_validate = temp[is.na(temp$validated), ]
    
    to_validate = to_validate[, c("person_id", "chrom", "start_pos", "end_pos",
        "ref_allele", "alt_allele", "hgnc", "consequence")]
    
    write.table(to_validate, file=path, quote=FALSE, sep="\t", row.names=FALSE)
}

main <- function() {
    
    validations = read.table(VALIDATION_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # get the genes with P values from the different test subsets
    genes = load_genes(ALL_PROBANDS_RESULTS, ALL_PROBANDS_RESULTS_WITHOUT_DIAGNOSED)
    
    genes$ddd_count = sapply(genes$hgnc, count_in_ddd, functional=FALSE)
    genes$ddd_count_without_validated = sapply(genes$hgnc, count_in_ddd, functional=FALSE, validated=validations)
    
    get_variants_to_validate(genes, validations,
        path="de_novos_for_validation_from_modified_indel_rate.tsv")
    
    
    # examine the de novos for each top-ranked gene one by one
    de_novos = get_raw_de_novos(DE_NOVOS_PATH)
    hgnc = "SLC6A1"
    examine_de_novos(hgnc, de_novos, validations)
}

main()
