#' Gene information.
#'
#' Table of all genes in the genome with HGNC symbol, chromosome and coding
#' sequence (CDS) length in base-pairs.
#'
#' @source url(http://www.ensembl.org/biomart)
#' @format A data frame with three variables: \code{hgnc}, \code{chrom},
#'     and \code{cds_length}
"gene_info"

#' Mutation rates for each mutation category for each gene.
#'
#' Table of all genes in the genome with refseq transcript ID, HGNC symbol,
#' length, then columns of mutation rates for different consequence types:
#'   a) log10-scaled mutuation rate for all consequence types
#'   b) log10-scaled mutuation rate for synonymous mutations
#'   c) log10-scaled mutuation rate for missense mutations
#'   e) log10-scaled mutuation rate for nonsense mutations
#'   f) log10-scaled mutuation rate for splice site mutations
#'   g) log10-scaled mutuation rate for frameshift mutations
#'
#' These mutation rates have been obtained from Nature Genetics 46:944-950
#' doi:10.1038/ng.3050 The only change to the table was to rename the "gene"
#' column to "hgnc".
#'
#' @source url(http://www.nature.com/ng/journal/v46/n9/extref/ng.3050-S2.xls)
#' @format A data frame with three variables: \code{transcript}, \code{hgnc},
#'         \code{bp}, \code{}, \code{syn}, \code{mis}, \code{non},
#'         \code{splice_site} and \code{frameshift}.
"gene_rates"
