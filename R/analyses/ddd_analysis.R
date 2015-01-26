# script to analyse enrichment of de novo mutations in genes in probands

library(mupit)

RATES_PATH = "/nfs/users/nfs_j/jm33/apps/de_novo_clustering/results/de_novo_gene_rates.ddd_4k.meta-analysis.without_chroms.txt"

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
    names(rates) = c("hgnc", "length", "mis", "non", "css", "syn", "frameshift")
    
    return(rates)
}

main <- function() {
    # # here's an example of how to use the functions in this script
    diagnosed_path = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx"
    # diagnosed = get_ddd_diagnosed(diagnosed_path)
    diagnosed = get_likely_diagnosed(diagnosed_path)
    # diagnosed = list(id=c(), sex=c())
    
    if (length(diagnosed$id) > 0) {
        prefix = "de_novos.ddd_4k.without_diagnosed"
    } else {
        prefix = "de_novos.ddd_4k.with_diagnosed"
    }
    
    rates = get_rates_dataset(RATES_PATH)
    
    # analyse the DDD only de novos
    trios = get_trio_counts(diagnosed)
    de_novos = get_de_novos(diagnosed)
    enriched = analyse_gene_enrichment(de_novos, trios,
        plot_path=paste("results/", prefix, ".ddd_only.manhattan.pdf", sep=""),
        rates=rates)
    
    # write the results, so that we can combine other analyses
    write.table(enriched, file=file.path("results",
        paste(prefix, ".ddd_only.enrichment_results.txt", sep="")), sep="\t",
        row.names=FALSE, quote=FALSE)
    
    # write the set of de novos for clustering analysis
    write.table(de_novos[, c("hgnc", "chrom", "start_pos", "consequence", "type")],
        file=paste(prefix, "ddd_only.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
    
    # analyse the DDD+ other cohorts de novos
    trios_meta = get_trio_counts(diagnosed, meta=TRUE)
    de_novos_meta = get_de_novos(diagnosed, meta=TRUE)
    enriched_meta = analyse_gene_enrichment(de_novos_meta, trios_meta,
        plot_path=paste("results/", prefix, ".meta-analysis.manhattan.pdf", sep=""),
        rates=rates)
    
    # write the results, so that we can combine other analyses
    write.table(enriched_meta, file=file.path("results",
        paste(prefix, ".meta-analysis.enrichment_results.txt", sep="")), sep="\t",
        row.names=FALSE, quote=FALSE)
    
    # write the set of de novos for clustering analysis
    write.table(de_novos_meta[, c("hgnc", "chrom", "start_pos", "consequence", "type")],
        file=paste(prefix, "meta-analysis.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
    
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
