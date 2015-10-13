# script to analyse enrichment of de novo mutations in genes in probands

library(argparse)
library(mupit)

parser = ArgumentParser()
parser$add_argument("--rates", help="Path to table of mutation rates.",
    default="/nfs/users/nfs_j/jm33/apps/denovonear/results/de_novo_gene_rates.ddd_4k.meta-analysis.txt")
parser$add_argument("--de-novos", help="Path to DDD de novo dataset.",
    default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-10-12.txt")
parser$add_argument("--validations", help="Path to validation results.",
    default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-10-12.txt")
parser$add_argument("--families", help="Path to families PED file.",
    default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt")
parser$add_argument("--trios", help="Path to file listing complete trios.",
    default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/trios.txt")
parser$add_argument("--ddg2p", help="Path to DDG2P file.",
    default="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter")
parser$add_argument("--diagnosed", help="Path to diagnosed probands file.")
parser$add_argument("--meta-analysis", default=FALSE, action="store_true",
    help="Whether to run meta-analysis that includes other published de novo datasets.")
parser$add_argument("--out-manhattan", help="Path to put PDF of manhattan plot.")
parser$add_argument("--out-probands-by-gene", help="Path to put json file of probands per gene.")
parser$add_argument("--out-enrichment", help="Path to put file of enrichment testing results.")
parser$add_argument("--out-clustering", help="Path to put file of enrichment testing results.")

#' defines the cohort sizes, used to get the overall population size
#'
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' @param meta true/false for whether to include meta-analysis populations
#'
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(families_path, trios_path, diagnosed_path, ddg2p_path, meta=FALSE) {
    
    families = read.table(families_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    families$is_proband = families$dad_id != "0" | families$mum_id != "0"
    
    # determine the trios with exome data available
    trios = read.table(trios_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    families$has_parental_data = families$individual_id %in% trios$proband_stable_id
    probands = families[families$is_proband & families$has_parental_data, ]
    
    # get the number of trios studied in our data for each sex
    sex = table(probands$sex)
    male = sex[["M"]]
    female = sex[["F"]] # female probands
    
    if (meta) {
        male = male + sum(publishedDeNovos::cohorts$unique_male)
        female = female + sum(publishedDeNovos::cohorts$unique_female)
    }
    
    # remove diagnosed patients, if maximising power
    if (!is.null(diagnosed_path)) {
        ddd_diagnosed = read.table(diagnosed_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        ddd_diagnosed = ddd_diagnosed[!duplicated(ddd_diagnosed[, c("person_id", "sex")]), ]
        
        ddg2p = load_ddg2p(ddg2p_path)
        external = publishedDeNovos::variants
        external_diagnosed = external[external$hgnc %in% ddg2p$gene[ddg2p$dominant] |
            (external$sex == "male" & external$hgnc %in% ddg2p$gene[ddg2p$hemizygous]), ]
        external_diagnosed = external_diagnosed[!duplicated(external_diagnosed[, c("person_id", "sex")]), ]
        
        # decrement for the diagnosed DDD individuals of each sex
        male = male - sum(ddd_diagnosed$sex %in% c("Male", "male", "M", "m"))
        female = female - sum(ddd_diagnosed$sex %in% c("Female", "female", "F", "f"))
        
        # decrement for the diagnosed external individuals of each sex
        male = male - sum(external_diagnosed$sex == "male", na.rm=TRUE)
        female = female - sum(external_diagnosed$sex == "female", na.rm=TRUE)
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
get_de_novos <- function(de_novos_path, validations_path, diagnosed_path, ddg2p_path, meta=FALSE) {
    
    diagnosed = NULL
    if (!is.null(diagnosed_path)) {
        diagnosed = read.table(diagnosed_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    }
    variants = get_ddd_de_novos(de_novos_path, diagnosed)
    
    validations = read.table(validations_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    variants = merge(variants, validations, by=c("person_id", "chrom",
        "start_pos", "end_pos", "ref_allele", "alt_allele", "hgnc",
        "consequence"), all.x=TRUE)
    
    # drop out the variants that failed to validate (i.e. were false positives,
    # or inherited)
    variants = variants[!variants$status %in% c("false_positive", "inherited"), ]
    variants$status = NULL
    
    if (meta) {
        external = publishedDeNovos::variants
        if (!is.null(diagnosed_path)) {
            ddg2p = load_ddg2p(ddg2p_path)
            external = external[!(external$hgnc %in% ddg2p$gene[ddg2p$dominant] |
                (external$sex == "male" & external$hgnc %in% ddg2p$gene[ddg2p$hemizygous])), ]
        }
        
        variants = rbind(variants, external)
    }
    
    return(variants)
}

get_rates_dataset <- function(rates_path) {
    rates = read.table(rates_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # convert from my column names to those used when estimating the gene
    # mutation rates given the cohort size
    rates$hgnc = rates$transcript_id
    rates$mis = rates$missense_rate
    rates$non = rates$nonsense_rate
    rates$splice_site = rates$splice_lof_rate
    rates$syn = rates$synonymous_rate
    rates$frameshift = rates$frameshift_rate
    
    rates = rates[, c("hgnc", "chrom", "length", "mis", "non", "splice_site", "syn", "frameshift")]
    
    return(rates)
}


#' run enrichment testing, and writes results as appropriate
#'
#' @param diagnosed list of sample IDs and sexes for diagnosed individuals
#' @param meta boolean for whether to include data for meta-analysis
run_tests <- function(de_novo_path, validations_path, families_path, trios_path,
    rates, diagnosed_path, ddg2p_path, meta, plot_path, json_path,
    enrichment_path, clustering_path) {
    
    # analyse the de novos
    trios = get_trio_counts(families_path, trios_path, diagnosed_path, ddg2p_path, meta)
    de_novos = get_de_novos(de_novo_path, validations_path, diagnosed_path, ddg2p_path, meta)
    enriched = mupit::analyse_gene_enrichment(de_novos, trios, plot_path=plot_path, rates=rates)
    
    # write the enrichment results to a table
    if (!is.null(enrichment_path)) {
        write.table(enriched, file=enrichment_path, sep="\t", row.names=FALSE, quote=FALSE)
    }
    
    # and write a list of probands with de novos per gene to a file. This is
    # for HPO similarity testing, so can only be used with DDD samples, since we
    # only have HPO phenotypes available for those individuals.
    if (!meta & !is.null(json_path)) { mupit::write_probands_by_gene(de_novos, json_path) }
    
    if (!is.null(clustering_path)) {
        # write the set of de novos for clustering analysis
        write.table(de_novos[, c("hgnc", "chrom", "start_pos", "consequence", "type")],
            file=clustering_path, sep="\t", row.names=FALSE, quote=FALSE)
    }
}

main <- function() {
    args = parser$parse_args()
    rates = get_rates_dataset(args$rates)
    
    run_tests(args$de_novos, args$validations, args$families, args$trios, rates,
        args$diagnosed, args$ddg2p, args$meta_analysis, args$plot_path,
        args$json_path, args$enrichment_path, args$clustering)
}


main()
