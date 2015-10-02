# script to analyse enrichment of de novo mutations in genes in probands

library(mupit)

PROBANDS_JSON_PREFIX = "/nfs/users/nfs_j/jm33/apps/mupit/data-raw/probands_by_gene"
RATES_PATH = "/nfs/users/nfs_j/jm33/apps/denovonear/results/de_novo_gene_rates.ddd_4k.meta-analysis.txt"
DE_NOVOS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-09-02.txt"
VALIDATIONS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-09-22.txt"
DIAGNOSED_PATH = "/lustre/scratch113/projects/ddd/users/jm33/ddd_4k.diagnosed.2015-10-02.txt"
DDG2P_PATH = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt"
TRIOS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/trios.txt"

#' defines the cohort sizes, used to get the overall population size
#'
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' @param meta true/false for whether to include meta-analysis populations
#'
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(families_path, trios_path, diagnosed_path, ddg2p_path, meta=FALSE) {
    
    families = read.table(families_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    families$is_proband = families$dad_id != "0" | families$mum_id != "0"
    
    # determine the trios with exome data available
    trios = read.table(trios_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    families$has_parental_data = families$individual_id %in% trios$proband_stable_id
    probands = families[families$is_proband & families$has_parental_data, ]
    
    # get the number of trios studied in our data for each sex
    sex = table(probands$sex)
    male = sex[["M"]]
    female = sex[["F"]] # female probands
    
    if (meta) {
        male = male + sum(publishedDeNovos::cohorts$unique_male)
        female = female + sum(publishedDeNovos::cohorts$unique_female)
    }
    
    # remove diagnosed patients, if maximising power
    if (!is.null(diagnosed_path)) {
        ddd_diagnosed = read.table(diagnosed_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        ddd_diagnosed = ddd_diagnosed[!duplicated(ddd_diagnosed[, c("person_id", "sex")]), ]
        
        ddg2p = load_ddg2p(ddg2p_path)
        external = publishedDeNovos::variants
        external_diagnosed = external[external$hgnc %in% ddg2p$gene[ddg2p$dominant] |
            (external$sex == "male" & external$hgnc %in% ddg2p$gene[ddg2p$hemizygous]), ]
        external_diagnosed = external_diagnosed[!duplicated(external_diagnosed[, c("person_id", "sex")]), ]
        
        # decrement for the diagnosed DDD individuals of each sex
        male = male - sum(ddd_diagnosed$sex %in% c("Male", "male", "M", "m"))
        female = female - sum(ddd_diagnosed$sex %in% c("Female", "female", "F", "f"))
        
        # decrement for the diagnosed external individuals of each sex
        male = male - sum(external_diagnosed$sex == "male", na.rm=TRUE)
        female = female - sum(external_diagnosed$sex == "female", na.rm=TRUE)
    }
    
    return(list(male=male, female=female))
}

#' combine datasets listing de novo mutations into a single data frame
#'
#' @param diagnosed list of IDs and sex for probands with diagnoses in the DDD
#' @param meta true/false for whether to include meta-analysis populations
#'
#' @return data frame containing HGNC, chrom, position, consequence, SNV or
#'    INDEL type, and study ID.
get_de_novos <- function(de_novos_path, validations_path, diagnosed_path, ddg2p_path, meta=FALSE) {
    
    diagnosed = NULL
    if (!is.null(diagnosed_path)) {
        diagnosed = read.table(diagnosed_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    }
    variants = get_ddd_de_novos(de_novos_path, diagnosed)
    
    validations = read.table(validations_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    variants = merge(variants, validations, by=c("person_id", "chrom",
        "start_pos", "end_pos", "ref_allele", "alt_allele", "hgnc",
        "consequence"), all.x=TRUE)
    
    # drop out the variants that failed to validate (i.e. were false positives,
    # or inherited)
    variants = variants[!variants$status %in% c("false_positive", "inherited"), ]
    variants$status = NULL
    
    if (meta) {
        external = publishedDeNovos::variants
        if (!is.null(diagnosed_path)) {
            ddg2p = load_ddg2p(ddg2p_path)
            external = external[!(external$hgnc %in% ddg2p$gene[ddg2p$dominant] |
                (external$sex == "male" & external$hgnc %in% ddg2p$gene[ddg2p$hemizygous])), ]
        }
        
        variants = rbind(variants, external)
    }
    
    return(variants)
}

get_rates_dataset <- function(rates_path) {
    rates = read.table(rates_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # convert from my column names to those used when estimating the gene
    # mutation rates given the cohort size
    rates$hgnc = rates$transcript_id
    rates$mis = rates$missense_rate
    rates$non = rates$nonsense_rate
    rates$splice_site = rates$splice_lof_rate
    rates$syn = rates$synonymous_rate
    rates$frameshift = rates$frameshift_rate
    
    rates = rates[, c("hgnc", "chrom", "length", "mis", "non", "splice_site", "syn", "frameshift")]
    
    return(rates)
}


#' run enrichment testing, and writes results as appropriate
#'
#' @param diagnosed list of sample IDs and sexes for diagnosed individuals
#' @param meta boolean for whether to include data for meta-analysis
run_tests <- function(de_novo_path, validations_path, families_path, trios_path, rates, diagnosed_path, ddg2p_path, meta) {
    prefix = "de_novos.ddd_4k.with_diagnosed"
    json_path = paste(PROBANDS_JSON_PREFIX, ".with_diagnosed.json", sep="")
    if (!is.null(diagnosed_path)) {
        prefix = "de_novos.ddd_4k.without_diagnosed"
        json_path = paste(PROBANDS_JSON_PREFIX, ".without_diagnosed.json", sep="")
    }
    
    # sometimes we run with DDD only samples, othertimes we want to include
    # the data from other cohorts for a meta-analaysis.
    mid_string = "ddd_only"
    if (meta) { mid_string = "meta-analysis" }
    
    # analyse the de novos
    trios = get_trio_counts(families_path, trios_path, diagnosed_path, ddg2p_path, meta)
    de_novos = get_de_novos(de_novo_path, validations_path, diagnosed_path, ddg2p_path, meta)
    enriched = mupit::analyse_gene_enrichment(de_novos, trios,
        plot_path=file.path("results", paste(prefix, mid_string, "manhattan",
        Sys.Date(), "pdf", sep=".")), rates=rates)
    
    # write the enrichment results to a table
    write.table(enriched, file=file.path("results",
        paste(prefix, mid_string, "enrichment_results", Sys.Date(), "txt",
        sep=".")), sep="\t", row.names=FALSE, quote=FALSE)
    
    # and write a list of probands with de novos per gene to a file. This is
    # for HPO similarity testing, so can only be used with DDD samples, since we
    # only have HPO phenotypes available for those individuals.
    if (!meta) { mupit::write_probands_by_gene(de_novos, json_path) }
    
    # write the set of de novos for clustering analysis
    write.table(de_novos[, c("hgnc", "chrom", "start_pos", "consequence", "type")],
        file=file.path("data-raw", paste(prefix, mid_string, "txt", sep=".")),
        sep="\t", row.names=FALSE, quote=FALSE)
}

main <- function() {
    rates = get_rates_dataset(RATES_PATH)
    
    run_tests(DE_NOVOS_PATH, VALIDATIONS_PATH, FAMILIES_PATH, TRIOS_PATH, rates, DIAGNOSED_PATH, DDG2P_PATH, meta=FALSE)
    run_tests(DE_NOVOS_PATH, VALIDATIONS_PATH, FAMILIES_PATH, TRIOS_PATH, rates, DIAGNOSED_PATH, DDG2P_PATH, meta=TRUE)
    
    DIAGNOSED_PATH = NULL
    run_tests(DE_NOVOS_PATH, VALIDATIONS_PATH, FAMILIES_PATH, TRIOS_PATH, rates, DIAGNOSED_PATH, DDG2P_PATH, meta=FALSE)
    run_tests(DE_NOVOS_PATH, VALIDATIONS_PATH, FAMILIES_PATH, TRIOS_PATH, rates, DIAGNOSED_PATH, DDG2P_PATH, meta=TRUE)
    
    # # Plot QQ plots for the meta-analysis de novos (requires statistics for
    # # all genes).
    # all_genes = analyse_gene_enrichment(de_novos_meta, trios, all_genes=TRUE)
    # qqman::qq(all_genes$p_func, main="Functional P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # qqman::qq(all_genes$p_lof, main="LoF P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # qqman::qq(all_genes$p_synonymous, main="Synomymous P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # dev.off()
}


main()
