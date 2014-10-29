# plotting functions for enrichment analyses

#' make plot labels for genes with fdr > threshold
#' 
#' @param enriched data frame containing HGNC symbols
#' @param p_values vector of p-values, sorted as per enriched data frame
#' @param num.tests number of tests performed (used for multiple correction).
#' @export
label_genes <- function(enriched, p_values, num.tests) {    
    fdr = p.adjust(p_values, method="BH", n=num.tests)
    fdr.thresh = 0.05
    
    # somethimes none of the genes have P values below the FDR threshold, return
    # from the function, rather than trying to continue
    if (!any(fdr < fdr.thresh)) {
        return()
    }
     
    label.thresh = max(p_values[which(fdr < fdr.thresh)])
    thresh.index = which(p_values < label.thresh)
    num.thresh = length(thresh.index)
    
    # sometimes we don't have any genes with P values more significant than the 
    # FDR threshold, simply return, rather than trying to plot (also since when
    # num.thresh is zero, we get some zero-length errors)
    if (num.thresh == 0) {
        return()
    }
    
    for (i in 1:num.thresh) {
        text(x=thresh.index[i], -log10(p_values[thresh.index[i]]), labels=enriched$hgnc[thresh.index[i]], pos=3, cex=0.5)
    }
}

#' make Manhattan plots for LOF and Func variants separately
#' 
#' @param num.tests number of tests performed (used by Bonferroni and FDR 
#'     correction).
#' @param enriched data frame containing HGNC symbols
#' @param length_p_vals vector of P values from testing using mutation rates 
#'     from length based estimates
#' @param daly_p_vals vector of P values from testing with mutation rates 
#'     provided by Mark Daly
#' @param title character string for the title of the plot
#' @param color_index vector of colors, for alternating colors per chromosome
#' @export
plot_values <- function(num.tests, enriched, length_p_vals, daly_p_vals, title, color_index) {
    # plot the results from de novos
    
    # set up the plot, starting with the length based P values
    plot(-log10(length_p_vals), col=color_index, pch=19, cex=0.75, 
        ylab="-log10(p)", xaxt="n", main=title, xlab="genome position",  
        ylim=c(0, max(-log10(length_p_vals), na.rm=TRUE) + 1))
    
    # plot line showing where Bonferroni correction would be
    abline(h = -log10(0.05 / num.tests), col="red", lty=2)
    
    # add the results from using the alternative mutation rates
    points(-log10(daly_p_vals), col=color_index, pch=1, cex=0.75)
    legend("topright", legend=c("Length-based rates", "Daly group rates"), 
        pch=c(19, 1), col="darkblue", cex=0.75)
    
    # label genes which are significant following FDR correction
    label_genes(enriched, length_p_vals, num.tests)
}

#' make Manhattan plots for LOF and Func variants separately
#' 
#' @param enriched data frame containing columns for chr, coord (position), and
#'         p values from testing for enrichment with loss-of-function and 
#'         functional consequence variants. Both of these have been tested 
#'         with two sets of mutation rate data, one derived from length 
#'         based rates, and the other from mutation rates provided by Mark 
#'         Daly.
#' @param num.tests number of tests performed (used by Bonferroni and FDR 
#'         correction).
#' @export
plot_enrichment_graphs <- function(enriched, num.tests) {
    
    Cairo::Cairo(file="temp.pdf", type="pdf", width=20, height=20, units = "cm")
    
    enriched = na.omit(enriched)
    
    # set up alternating colors for successive chromosomes
    enriched$chrom = as.character(enriched$chrom)
    enriched$chrom[enriched$chrom == "X"] = "23"
    enriched = enriched[order(as.numeric(as.character(enriched$chrom)), as.numeric(as.character(enriched$min_pos))), ]
    color_index = rep("lightblue3", nrow(enriched))
    odd.chr = c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19", "21", "23")
    color_index[enriched$chrom %in% odd.chr] = "darkblue"
    
    # plot the results from loss-of-function de novos
    plot_values(num.tests, enriched, enriched$p.lof.length, enriched$p.lof.daly, "Loss-of-Function DNMs", color_index)
    plot_values(num.tests, enriched, enriched$p.func.length, enriched$p.func.daly, "Functional DNMs", color_index)
    
    dev.off()
}
