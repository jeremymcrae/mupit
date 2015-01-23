# script to analyse enrichment of de novo mutations in genes in probands

library(mupit)

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

main <- function() {
    # # here's an example of how to use the functions in this script
    diagnosed_path = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx"
    # diagnosed = get_ddd_diagnosed(diagnosed_path)
    diagnosed = get_likely_diagnosed(diagnosed_path)
    
    no_diagnosed = list(id=c(), sex=c())
    
    # analyse the DDD only de novos
    trios = get_trio_counts(diagnosed)
    de_novos = get_de_novos(diagnosed)
    enriched = analyse_gene_enrichment(de_novos, trios, "results/de_novos.ddd_4k.ddd_only.manhattan.pdf")
    
    # write the results, so that we can combine other analyses
    write.table(enriched, file=file.path("results",
        "de_novos.ddd_4k.without_diagnosed.ddd_only.enrichment_results.txt"), sep="\t",
        row.names=FALSE, quote=FALSE)
    
    # write the set of de novos for clustering analysis
    write.table(de_novos[, c("hgnc", "chrom", "start_pos", "consequence", "type")],
        file="de_novos.ddd_4k.ddd_only.txt", sep="\t", row.names=FALSE, quote=FALSE)
    
    # analyse the DDD+ other cohorts de novos
    trios_meta = get_trio_counts(diagnosed, meta=TRUE)
    de_novos_meta = get_de_novos(diagnosed, meta=TRUE)
    enriched_meta = analyse_gene_enrichment(de_novos_meta, trios_meta, "results/de_novos.ddd_4k.meta-analysis.manhattan.pdf")
    
    # write the results, so that we can combine other analyses
    write.table(enriched_meta, file=file.path("results",
        "de_novos.ddd_4k.without_diagnosed.meta-analysis.enrichment_results.txt"), sep="\t",
        row.names=FALSE, quote=FALSE)
    
    # write the set of de novos for clustering analysis
    write.table(de_novos_meta[, c("hgnc", "chrom", "start_pos", "consequence", "type")],
        file="de_novos.ddd_4k.meta-analysis.txt", sep="\t", row.names=FALSE, quote=FALSE)
    
    # # PLot QQ plots for the meta-analysis de novos (requires statistics for
    # # all genes).
    # all_genes = analyse_gene_enrichment(de_novos_meta, trios, all_genes=TRUE)
    # qqman::qq(all_genes$p_func, main="Functional P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # qqman::qq(all_genes$p_lof, main="LoF P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # qqman::qq(all_genes$p_synonymous, main="Synomymous P values", las=1, cex.lab=1.4, cex.axis=1.4)
    # dev.off()
    
    head(enriched[order(enriched$p_func), ])
}


main()
