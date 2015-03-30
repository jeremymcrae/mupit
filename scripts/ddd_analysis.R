# script to analyse enrichment of de novo mutations in genes in probands

library(mupit)

PROBANDS_JSON_PREFIX = "/nfs/users/nfs_j/jm33/apps/mupit/data-raw/probands_by_gene"
RATES_PATH = "/nfs/users/nfs_j/jm33/apps/de_novo_clustering/results/de_novo_gene_rates.ddd_4k.meta-analysis.txt"

#' defines the cohort sizes, used to get the overall population size
#'
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' @param meta true/false for whether to include meta-analysis populations
#'
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(diagnosed, meta=FALSE) {
    
    # number of trios studied in our data
    male = 2408 # male probands
    female = 1887 # female probands
    
    # remove diagnosed patients, if maximising power
    male = male - sum(diagnosed$sex %in% c("Male", "male", "M", "m"))
    female = female - sum(diagnosed$sex %in% c("Female", "female", "F", "f"))
    
    if (meta) {
        male = male + sum(cohorts$unique_male)
        female = female + sum(cohorts$unique_female)
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
get_de_novos <- function(diagnosed, meta=FALSE) {
    
    variants = get_ddd_de_novos(diagnosed)
    
    if (meta) {
        variants = rbind(variants, published_de_novos)
    }
    
    return(variants)
}

get_rates_dataset <- function(rates_path) {
    rates = read.table(rates_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    names(rates) = c("hgnc", "chrom", "length", "mis", "non", "css", "splice_region", "syn", "frameshift")
    
    return(rates)
}


#' run enrichment testing, and writes results as appropriate
#'
#' @param diagnosed list of sample IDs and sexes for diagnosed individuals
#' @param meta boolean for whether to include data for meta-analysis
run_tests <- function(rates, diagnosed, meta) {
    prefix = "de_novos.ddd_4k.with_diagnosed"
    json_path = paste(PROBANDS_JSON_PREFIX, ".with_diagnosed.json", sep="")
    if (length(diagnosed$id) > 0) {
        prefix = "de_novos.ddd_4k.without_diagnosed"
        json_path = paste(PROBANDS_JSON_PREFIX, ".without_diagnosed.json", sep="")
    }
    
    # sometimes we run with DDD only samples, othertimes we want to include
    # the data from other cohorts for a meta-analaysis.
    mid_string = "ddd_only"
    if (meta) { mid_string = "meta-analysis" }
    
    # analyse the de novos
    trios = get_trio_counts(diagnosed, meta)
    de_novos = get_de_novos(diagnosed, meta)
    enriched = mupit::analyse_gene_enrichment(de_novos, trios,
        plot_path=file.path("results", paste(prefix, mid_string, "manhattan.pdf", sep=".")),
        rates=rates)
    
    # write the enrichment results to a table
    write.table(enriched, file=file.path("results",
        paste(prefix, mid_string, "enrichment_results.txt", sep=".")), sep="\t",
        row.names=FALSE, quote=FALSE)
    
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
    
    # here's an example of how to use the functions in this script
    diagnosed_path = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx"
    diagnosed = get_likely_diagnosed(diagnosed_path)
    run_tests(rates, diagnosed, meta=FALSE)
    run_tests(rates, diagnosed, meta=TRUE)
    
    diagnosed = NULL
    run_tests(rates, diagnosed, meta=FALSE)
    run_tests(rates, diagnosed, meta=TRUE)
    
    # # Plot QQ plots for the meta-analysis de novos (requires statistics for
    # # all genes).
    # all_genes = analyse_gene_enrichment(de_novos_meta, trios, all_genes=TRUE)
    # qqman::qq(all_genes$p_func, main="Functional P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # qqman::qq(all_genes$p_lof, main="LoF P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # qqman::qq(all_genes$p_synonymous, main="Synomymous P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # dev.off()
}


main()
