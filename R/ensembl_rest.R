# functions to extract data from the ensembl REST API


# consequence list, as sorted at http://www.ensembl.org/info/genome/variation/predicted_data.html
consequences = c("transcript_ablation", "splice_donor_variant", "splice_acceptor_variant", "stop_gained", "frameshift_variant", "stop_lost", "initiator_codon_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "regulatory_region_variant", "feature_elongation", "feature_truncation", "intergenic_variant")
severity = data.frame(consequence=consequences, rank=seq(1:length(consequences)),
    stringsAsFactors=FALSE)

#' find the VEP consequence for a variant
#' 
#' @param variant data frame or list for a variant, containing columns named
#'     "chrom", "start_pos", "end_pos", and "alt_allele" code for a single variant
#' @param build genome build to find consequences on
#' @param verbose flag indicating whether to print variants as they are checked
#' 
#' @export
#' @return a character string containing the most severe consequence, as per VEP
#'     annotation formats.
#' 
#' @examples
#' get_vep_consequence(data.frame(chrom=c("1"), 
#'     start_pos=c("1000000"), end_pos=c("1000000"), alt_allele=c("A"), 
#'     ref_allele=c("G")))
#' get_vep_consequence(list(chrom="1", start_pos="1000000", 
#'     end_pos="1000000", alt_allele="A", ref_allele="G"))
get_vep_consequence <- function(variant, build="grch37", verbose=FALSE) {
    
    # only tolerate the grch37 and grch38 genome builds, since they are the only
    # genome builds supported by the Ensembl REST API
    allowed_builds = c("grch37", "grch38")
    stopifnot( tolower(build) %in% allowed_builds )
    
    # get the correct allele for deletions for the Ensembl REST API requests
    if (variant[["alt_allele"]] == "") {variant[["alt_allele"]] = "-"}
    
    base_url = "rest.ensembl.org/vep/human/region/"
    if (build == "grch37") {
        base_url = paste("grch37", base_url, sep = ".")
    }
    
    # define parts of the URL
    url_start = paste("http://", base_url, variant[["chrom"]], ":", 
        variant[["start_pos"]], ":", variant[["end_pos"]], "/", sep="")
    content = "?content-type=application/json"
    
    # define the URLs for the ref and alt alleles
    url = paste(url_start, variant[["alt_allele"]], content, sep = "")
    ref_url = paste(url_start, variant[["ref_allele"]], content, sep = "")
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], " ", 
            variant[["alt_allele"]], "    ", url, sep=""))
    }
    
    json = try(rjson::fromJSON(file=url), silent=TRUE)
    if (class(json) == "try-error") {json = rjson::fromJSON(file=ref_url)}
    
    transcript = find_most_severe_transcript(json)
    
    consequence = transcript$consequence_terms
    temp_severity = min(severity$rank[severity$consequence %in% consequence])
    consequence = severity$consequence[severity$rank == temp_severity]
    
    value = list(consequence=consequence, gene=transcript$gene_symbol)
    
    return(value)
}

#' find the VEP consequence for a variant
#' 
#' @param ensembl_json json data for variant from Ensembl
#' 
#' @export
#' @return a character string containing the most severe consequence, as per VEP
#'     annotation formats.
find_most_severe_transcript <- function(ensembl_json) {
    
    bad_transcripts = c("lincRNA", "nonsense_mediated_decay", "transcribed_unprocessed_pseudogene")
    
    best_transcript = NA
    best_severity = NA
    
    for (transcript in ensembl_json[[1]]$transcript_consequences) {
        
        # don't bother to check some transcript types
        if (transcript$biotype %in% bad_transcripts) { next }
        
        # get consequence and severity rank in the current transcript
        consequence = transcript$consequence_terms
        temp_severity = min(severity$rank[severity$consequence %in% consequence])
        consequence = severity$consequence[severity$rank == temp_severity]
        
        # check if this is the most severe consequence
        if (is.na(best_severity) | temp_severity < best_severity) {
            best_severity = temp_severity
            best_transcript = transcript
        }
    }
    
    return(best_transcript)
}

#' find the hgnc symbol overlapping a variant position
#' 
#' @param variant data frame or list for a variant, containing columns named
#'     "chrom", "start_pos", and "end_pos" for a single variant
#' @param build genome build to find consequences on
#' @param verbose flag indicating whether to print variants as they are checked
#' 
#' @export
#' @return a character string containing the HGNC symbol.
#' 
#' @examples
#' get_gene_id_for_variant(data.frame(chrom=c("1"), start_pos=c("1000000"),
#'     end_pos=c("1000000")))
#' get_gene_id_for_variant(list(chrom="1", start_pos="1000000", 
#'     end_pos="1000000"))
get_gene_id_for_variant <- function(variant, build="grch37", verbose=FALSE) {
    
    # only tolerate the grch37 and grch38 genome builds, since they are the only
    # genome builds supported by the Ensembl REST API
    allowed_builds = c("grch37", "grch38")
    stopifnot( tolower(build) %in% allowed_builds )
    
    base_url = "rest.ensembl.org/overlap/region/human/"
    if (build == "grch37") {
        base_url = paste("grch37", base_url, sep = ".")
    }
    
    # make sure the end position is suitable for the Ensembl REST API request
    if (as.numeric(variant[["end_pos"]]) < as.numeric(variant[["start_pos"]])) {
        variant[["end_pos"]] = variant[["start_pos"]]
    }
    
    # define parts of the URL
    url_start = paste("http://", base_url, variant[["chrom"]], ":", 
        variant[["start_pos"]], ":", variant[["end_pos"]], "/", sep="")
    content = "?feature=gene;content-type=application/json"
    url = paste(url_start, content, sep = "")
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], 
            "    ", url, sep=""))
    }
    
    json = rjson::fromJSON(file=url)
    
    # return blank string for variants not in genes
    if (length(json) > 0) {
        return(json[[1]]$external_name)
    }
    
    return("")
}

#' find genomic sequence within a region
#' 
#' @param variant data frame or list for a variant, containing columns named
#'     "chrom", "start_pos", and "end_pos" for a single variant
#' @param build genome build to find consequences on
#' @param verbose flag indicating whether to print variants as they are checked
#' 
#' @export
#' @return a character string containing the HGNC symbol.
#' 
#' @examples
#' get_gene_id_for_variant(data.frame(chrom=c("1"), start_pos=c("1000000"),
#'     end_pos=c("1000000")))
#' get_gene_id_for_variant(list(chrom="1", start_pos="1000000", 
#'     end_pos="1000000"))
get_sequence_in_region <- function(variant, build="grch37", verbose=FALSE) {
    
    # only tolerate the grch37 and grch38 genome builds, since they are the only
    # genome builds supported by the Ensembl REST API
    allowed_builds = c("grch37", "grch38")
    stopifnot( tolower(build) %in% allowed_builds )
    
    base_url = "rest.ensembl.org/sequence/region/human/"
    if (build == "grch37") {
        base_url = paste("grch37", base_url, sep = ".")
    }
    
    # make sure the end position is suitable for the Ensembl REST API request
    if (as.numeric(variant[["end_pos"]]) < as.numeric(variant[["start_pos"]])) {
        variant[["end_pos"]] = variant[["start_pos"]]
    }
    
    # define the URL
    url = paste("http://", base_url, variant[["chrom"]], ":", 
        variant[["start_pos"]], ":", variant[["end_pos"]], 
        ":1?content-type=application/json", sep="")
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], 
            "    ", url, sep=""))
    }
    
    json = rjson::fromJSON(file=url)
    
    # return blank string for variants not in genes
    if (length(json) > 0) {
        return(json$seq)
    }
    
    return("")
}
