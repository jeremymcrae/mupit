#' De novos from published studies.
#'
#' De novo variants found within participants in published studies. The studies 
#' are:
#' 
#' Rauch et al. (2012) Lancet 380:1674-1682. Supplementary tables 2 and 3
#' doi: 10.1016/S0140-6736(12)61480-9
#' 
#' De Ligt et al. (2012) N Engl J Med 367:1921-1929. Supplementary table 3
#' doi: 10.1056/NEJMoa1206524
#'
#' Gilissen et al. (2014) Nature 511:344-347. Supplementary table 8
#' doi: 10.1038/nature13394
#'
#' American Journal of Human Genetics (2014) 95:360-370. Supplementary table 1
#' doi: 10.1016/j.ajhg.2014.08.013
#'
#' Sanders et al. (2012) Nature 485:237-241. Supplementary table 2
#' doi: 10.1038/nature10945
#'
#' O'Roak et al. (2012) Nature 485:246-250. Supplementary table 3
#' doi: 10.1038/nature10989
#'
#' Iossifov et al. (2012) Neuron 74:285-299. Supplementary tables 1 and 2
#' doi: 10.1016/j.neuron.2012.04.009
#'
#' Iossifov et al. (2014) Nature 498:216-221. Supplementary table 2
#' doi: 10.1038/nature13908
#'
#' De Rubeis et al. (2013) Nature 515:209-215. Supplementary table 3
#' doi: 10.1038/nature13772
#'
#' Fromer et al. (2014) Nature 506:179-184. Supplementary table 1
#' doi: 10.1038/nature12929
#'
#' Zaidi et al. (2013) Nature 498:220-223. Supplementary table 4
#' doi: 10.1038/nature12141
#'
#' @source url("http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0140673612614809/1-s2.0-S0140673612614809-mmc1.pdf/271074/FULL/S0140673612614809/55b26043f4a279334b3a5ec00b9faf4b/mmc1.pdf")
#' @source url("http://www.nejm.org/doi/suppl/10.1056/NEJMoa1206524/suppl_file/nejmoa126524_appendix.pdf")
#' @source url("http://www.nature.com/nature/journal/v511/n7509/extref/nature13394-s1.pdf")
#' @source url("http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0002929714003838/1-s2.0-S0002929714003838-mmc2.xlsx/276895/FULL/S0002929714003838/bf21945d72e3297fc44969dc0296f4f1/mmc2.xlsx")
#' @source url("http://www.nature.com/nature/journal/v485/n7397/extref/nature10945-s3.xls")
#' @source url("http://www.nature.com/nature/journal/v485/n7397/extref/nature10989-s2.xls")
#' @source url("http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0896627312003406/1-s2.0-S0896627312003406-mmc2.xlsx/272195/FULL/S0896627312003406/26c5ba3b72a2410ef43fec52a40f35e6/mmc2.xlsx")
#' @source url("http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0896627312003406/1-s2.0-S0896627312003406-mmc4.xlsx/272195/FULL/S0896627312003406/6caa42b35609c2ed5910b5381ddd5335/mmc4.xlsx")
#' @source url("http://www.nature.com/nature/journal/v515/n7526/extref/nature13908-s2.zip")
#' @source url("http://www.nature.com/nature/journal/v515/n7526/extref/nature13772-s4.xlsx")
#' @source url("http://www.nature.com/nature/journal/v506/n7487/extref/nature12929-s2.xlsx")
#' @source url("http://www.nature.com/nature/journal/v498/n7453/extref/nature12141-s1.pdf")
#' @format A data frame with twelve variables: c(\code{person_id}, \code{chrom},
#'     \code{start_pos}, \code{end_pos}, \code{ref_allele}, \code{alt_allele}, 
#'     \code{hgnc}, \code{consequence}, \code{study_code}, 
#'     \code{publication_doi}, \code{study_phenotype} and \code{type}
"published_de_novos"

#' Gene information.
#' 
#' Table of all genes in the genome with HGNC symbol, chromosome and coding
#' sequence (CDS) length in base-pairs.
#' 
#' @source url(http://www.ensembl.org/biomart)
#' @format A data frame with three variables: \code{hgnc}, \code{chrom},
#'     and \code{cds_length}
"gene_info"

#' Cohort information.
#' 
#' Table of sex-specific numbers of affected probands reported on in published
#' exome sequencing of children with developmental disorders (autism, congenital
#' heart disease, intellectual disability, epilepsy, and schizophrenia). 
#' 
#' This provides the total population of trio-based probands investigated within
#' these studies, that is, the population for whom a de novo mutation might have
#' been identified. This determines the number of mutations we expect to find
#' within the total population.
#' 
#' The table of sex-counts has been carefully extracted from the publications, 
#' accounting for overlap in samples between different publications.
#' 
#' @format A data frame with eight variables: \code{study_code}, \code{year}, 
#'     \code{study_phenotype}, \code{publication_doi}, \code{unique_male}, 
#'     \code{unique_female}, \code{subsets}, and \code{comment}.
"cohorts"