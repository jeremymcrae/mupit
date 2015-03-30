# plotting functions for enrichment analyses

#' make plot labels for genes with fdr > threshold
#'
#' @param enriched data frame containing HGNC symbols
#' @param p_values vector of p-values, sorted as per enriched data frame
#' @param num_tests number of tests performed (used for multiple correction).
#' @export
label_genes <- function(enriched, p_values, num_tests) {
    fdr = p.adjust(p_values, method="BH", n=num_tests)
    fdr.thresh = 0.05
    
    # somethimes none of the genes have P values below the FDR threshold, return
    # from the function, rather than trying to continue
    if (!any(fdr < fdr.thresh)) {
        return()
    }
    
    thresh.index = p_values < max(p_values[fdr < fdr.thresh])
    
    # sometimes we don't have any genes with P values more significant than the
    # FDR threshold, simply return, rather than trying to plot (also since when
    # num.thresh is zero, we get some zero-length errors)
    if (sum(thresh.index) == 0) {
        return()
    }
    
    text(x=which(thresh.index), y=-log10(p_values[thresh.index]),
        labels=enriched$hgnc[thresh.index], pos=3, cex=0.6)
}

#' make Manhattan plots for LoF, func and silent variants separately
#'
#' @param num_tests number of tests performed (used by Bonferroni and FDR
#'     correction).
#' @param enriched data frame containing HGNC symbols
#' @param p_values vector of P values from enrichment testing
#' @param title character string for the title of the plot
#' @export
plot_values <- function(num_tests, enriched, p_values, title) {
    # plot the results from de novos
    
    # set up alternating colors for successive chromosomes
    color_index = rep("lightblue3", nrow(enriched))
    odd.chr = c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19", "21", "23")
    color_index[enriched$chrom %in% odd.chr] = "darkblue"
    
    # set up the plot, starting with the length based P values
    plot(-log10(p_values), col=color_index, pch=19, cex=0.75,
        ylab=expression(-log[10](italic(P))), xaxt="n", main=title,
        xlab="genome position", ylim=c(0, max(-log10(p_values), na.rm=TRUE) + 1),
        las=1, cex.axis=1.4, cex.lab=1.4)
    
    # plot line showing where Bonferroni correction would be
    abline(h = -log10(0.05 / num_tests), col="red", lty=2)
    
    # label genes which are significant following FDR correction
    label_genes(enriched, p_values, num_tests)
}

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
    
    Cairo::Cairo(file=path, type="pdf", width=15, height=15, units = "cm")
    
    enriched = na.omit(enriched)
    
    # sort the dataset by chromosome
    enriched$chrom = as.character(enriched$chrom)
    enriched$chrom[enriched$chrom == "X"] = "23"
    enriched = enriched[order(as.numeric(enriched$chrom), as.numeric(as.character(enriched$min_pos))), ]
    
    # plot the results from loss-of-function and functional de novos
    plot_values(num_tests, enriched, enriched$p_lof, "Loss-of-Function DNMs")
    plot_values(num_tests, enriched, enriched$p_func, "Functional DNMs")
    
    dev.off()
}
