# plotting functions for enrichment analyses

#' make Manhattan plots for LOF and Func variants separately
#'
#' @param enriched data frame containing columns for chr, coord (position), and
#'         p values from testing for enrichment with loss-of-function and
#'         functional consequence variants.
#' @param num_tests number of tests performed (used by Bonferroni and FDR
#'         correction).
#' @param path path to save the plot to
#' @export
plot_enrichment_graphs <- function(enriched, num_tests, path) {
    
    Cairo::Cairo(file=path, type="pdf", width=20, height=15, units = "cm")
    
    colors = c("lightblue3", "darkblue")
    enriched = na.omit(enriched)
    
    # sort the dataset by chromosome
    enriched$chrom = as.character(enriched$chrom)
    enriched$chrom[enriched$chrom == "X"] = "23"
    enriched$chrom = as.numeric(enriched$chrom)
    
    # get the p-value to a FDR threshold of 0.05, below which
    # we will plot gene labels
    p_values = enriched$enrichment_p_value
    fdr = p.adjust(p_values, method="BH", n=num_tests)
    threshold = min(p_values[fdr >= 0.05])
    
    qqman::manhattan(enriched,
        chr="chrom", bp="min_pos", snp="hgnc", p="enrichment_p_value",
        chrlabs=c(1:22, "X"), annotatePval=threshold, annotateTop=FALSE,
        suggestiveline=FALSE, genomewideline=-log10(0.05/num_tests),
        col=colors, las=2)
    
    dev.off()
}
