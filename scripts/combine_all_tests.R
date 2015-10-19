# script to combine different test statistics together into a single table

library(argparse)
library(mupit)

get_options <- function() {
    parser = ArgumentParser()
    parser$add_argument("--ddg2p", help="Path to DDG2P file.")
    parser$add_argument("--ddd-enrichment",
        help="Path to results from enrichment testing of DDD de novos.")
    parser$add_argument("--meta-enrichment",
        help="Path to results from enrichment testing of meta-analysis de novos.")
    parser$add_argument("--ddd-clustering",
        help="Path to results from proximity clustering testing of DDD de novos.")
    parser$add_argument("--meta-clustering",
        help="Path to results from proximity clustering testing of meta-analysis de novos.")
    parser$add_argument("--ddd-phenotype",
        help="Path to file of DDD HPO similarity testing results.")
    parser$add_argument("--output",
        help="Path to write output file for combined testing results.")
    
    args = parser$parse_args()
    
    return(args)
}

include_ddg2p_status <- function(merged, ddg2p_path) {
    # load the DDG2P, so that we can annotated the genes with their DDG2P status
    ddg2p = load_ddg2p(ddg2p_path)
    
    # annotate each column with DDG2P status
    merged$in_ddg2p = merged$hgnc %in% ddg2p$gene
    merged$in_ddg2p_dominant = merged$hgnc %in% ddg2p$gene[ddg2p$dominant]
    
    return(merged)
}

main <- function() {
    args = get_options()
    
    merged = combine_tests(args$meta_clustering, args$meta_enrichment,
        args$ddd_clustering, args$ddd_enrichment, args$ddd_phenotype)
    merged = include_ddg2p_status(merged, args$ddg2p)
    write.table(merged, file=args$output, row.names=F, quote=F, sep="\t")
}

main()
