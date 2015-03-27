#' Write a json-formatted list of probands per gene
#'
#' We want a correctly-formatted JSON list of probands per gene, for analysis of
#' phenotype similarity between probands, see
#' \url{https://github.com/jeremymcrae/hpo_similarity}.
#'
#' @param de_novos dataframe of de-novo mutations per proband
#' @param path path to write the json data to
#' @export
write_probands_by_gene <- function(de_novos, path) {
    proband_genes = sort(unique(de_novos$hgnc))
    probands_json = sapply(proband_genes, function(x) de_novos$person_id[de_novos$hgnc == x])
    probands_json = jsonlite::toJSON(probands_json, pretty=TRUE)
    
    # write a json file of probands by gene, for HPO similarity testing
    fileConn = file(path)
    writeLines(probands_json, fileConn)
    close(fileConn)
}
