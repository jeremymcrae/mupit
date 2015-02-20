# functions to extract data from the ensembl REST API


# consequence list, as sorted at http://www.ensembl.org/info/genome/variation/predicted_data.html
consequences = c("transcript_ablation", "splice_donor_variant", "splice_acceptor_variant", "stop_gained", "frameshift_variant", "stop_lost", "initiator_codon_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "regulatory_region_variant", "feature_elongation", "feature_truncation", "intergenic_variant")
severity = data.frame(consequence=consequences, rank=seq(1:length(consequences)),
    stringsAsFactors=FALSE)

timeEnv = new.env()
assign("initial_time", Sys.time(), envir=timeEnv)

#' make a URL request to the Ensembl service
#'
#' @param ext extension for URL
#' @param build genome build for REST API
#' @param headers headers for the url
#' @param tries number of attempts that have been mader to access the URL
#'
#' @export
#' @return a character string, typically json encoded
#'
#' @examples
#' request_from_ensembl("/vep/human/region/1:205901016:205901016/A")
request_from_ensembl <- function(ext, headers="", build="grch37", tries=0) {
    
    base_url = "rest.ensembl.org"
    # strip any possible whitespace from the url, since imported datasets
    # occasionally contain whitespace in the values used to construct the URL.
    ext = gsub(" ", "", ext)
    
    # only tolerate the grch37 and grch38 genome builds, since they are the only
    # genome builds supported by the Ensembl REST API
    allowed_builds = c("grch37", "grch38")
    stopifnot( tolower(build) %in% allowed_builds )
    if (build == "grch37") {
        base_url = paste("grch37", base_url, sep = ".")
    }
    
    content = paste("?", paste(headers, "content-type=application/json", sep=";"), sep="")
    url = paste("http://", base_url, ext, content, sep="")
    
    # check that we are not requesting urls fater than that allowed by Ensembl,
    # sleep until the period expires
    current_time = Sys.time()
    diff = 0.067 - as.numeric(current_time - get("initial_time", current_time, envir=timeEnv))
    if (diff > 0) { Sys.sleep(diff) }
    
    # set the previous request time to that of the current request
    assign("initial_time", current_time, envir=timeEnv)
    
    # cut out after making 5 attempts to access the url
    tries = tries + 1
    stopifnot(tries <= 5)
    
    request = httr::GET(url)
    
    # handle the possible http request return status codes, such as when the
    # server is unavailable, when we have made too many requests, or requested
    # an impossible URL.
    if (request$status_code == 503) { # server down
        Sys.sleep(30)
        return(request_from_ensembl(url, tries))
    } else if (request$status_code == 429) { # too frequent requests
        reset_time = as.numeric(request$headers$`x-ratelimit-reset`)
        Sys.sleep(reset_time)
        return(request_from_ensembl(url, tries))
    } else if (request$status_code == 400) { # bad url
        msg = paste("bad url request: ", url, sep = "")
        stop(msg)
    } else if (request$status_code != 200) { # other status errors
        msg = paste("unknown status: ", request$status$code, " at: ", url, sep = "")
        stop(msg)
    }
    
    return(intToUtf8(request$content))
}

#' find the VEP consequence for a variant
#'
#' @param variant data frame or list for a variant, containing columns named
#'     "chrom", "start_pos", "end_pos", and "alt_allele" code for a single variant
#' @param include_hgvs whether to also include the HGVS consequence strings for
#'     cDNA and protein.
#' @param include_deleterious_score whether to also include SIFT and polyphen predictions
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
get_vep_consequence <- function(variant, include_hgvs=FALSE, include_deleterious_score=FALSE, build="grch37", verbose=FALSE) {
    
    # get the correct allele for deletions for the Ensembl REST API requests
    if (variant[["alt_allele"]] == "") {variant[["alt_allele"]] = "-"}
    
    # define parts of the URL
    ext = "/vep/human/region/"
    ext = paste(ext, variant[["chrom"]], ":", variant[["start_pos"]], ":", variant[["end_pos"]], "/", sep="")
    
    # define the URLs for the ref and alt alleles
    alt_ext = paste(ext, variant[["ref_allele"]], sep = "")
    ext = paste(ext, variant[["alt_allele"]], sep = "")
    
    headers = ""
    if (include_hgvs) {
        headers = "hgvs=1"
    }
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]], " ",
            variant[["alt_allele"]], "    ", ext, sep=""))
    }
    
    json = try(request_from_ensembl(ext, headers=headers, build))
    if (class(json) == "try-error") {json = request_from_ensembl(alt_ext, headers=headers, build)}
    json = rjson::fromJSON(json)
    
    transcript = find_most_severe_transcript(json, include_hgvs)
    
    consequence = transcript$consequence_terms
    temp_severity = min(severity$rank[severity$consequence %in% consequence])
    consequence = severity$consequence[severity$rank == temp_severity]
    
    value = list(consequence=consequence, gene=transcript$gene_symbol)
    if (include_hgvs) {
        value$hgvsc = transcript$hgvsc
        value$hgvsp = transcript$hgvsp
    }
    
    if (include_deleterious_score) {
        if ("polyphen_prediction" %in% names(transcript)) {
            value$polyphen_prediction = transcript$polyphen_prediction
            value$polyphen_score = transcript$polyphen_score
            value$sift_prediction = transcript$sift_prediction
            value$sift_score = transcript$sift_score
        } else {
            value$polyphen_prediction = NA
            value$polyphen_score = NA
            value$sift_prediction = NA
            value$sift_score = NA
        }
    }
    
    return(value)
}

#' find the most severe transcript from Ensembl data
#'
#' @param ensembl_json json data for variant from Ensembl
#' @param include_hgvs boolean showing if we will be using HGVS predictions from
#'     the VEP consequence predictions. If true, we sort the trasncripts by their
#'     length, so that we consistently return the same transcript between
#'     different variants from the same gene. Otherwise we just return the first
#'     most severe transcript.
#' @param exclude_bad boolean to show whether we want to exclude nonsense
#'     mediated decay transcripts and so forth. We only supply a false value for
#'     recursive calling, when we have failed to find any transcript during
#'     normal usage.
#'
#' @export
#' @return a character string containing the most severe consequence, as per VEP
#'     annotation formats.
find_most_severe_transcript <- function(ensembl_json, include_hgvs=FALSE, exclude_bad=TRUE) {
    
    bad_transcripts = c("lincRNA", "nonsense_mediated_decay", "pseudogene",
        "transcribed_unprocessed_pseudogene")
    
    best_transcript = NA
    best_severity = NA
    symbol = NA
    
    # if we are dealing with an intergenic variant, the variant won't have any
    # transcript consequences, and so the function will enter an infinite
    # recursion. Instead return the json entry (appending the correct
    # annotations for the higher function)
    if (!("transcript_consequences" %in% names(ensembl_json[[1]]))) {
        nontranscript = ensembl_json[[1]]
        nontranscript$consequence_terms = nontranscript$most_severe_consequence
        nontranscript$gene_symbol = ""
        return(nontranscript)
    }
    
    transcripts = ensembl_json[[1]]$transcript_consequences
    # sort the transcripts by their protein lengths
    if (include_hgvs) {
        protein_lengths = sapply(transcripts, function(x) get_protein_length(x$transcript_id))
        protein_order = sort(protein_lengths, decreasing=TRUE, index.return=TRUE)
        transcripts = lapply(protein_order$ix, function(x) transcripts[[x]])
    }
    
    for (transcript in transcripts) {
        
        # don't bother to check some transcript types, depending on whether we
        if (exclude_bad & transcript$biotype %in% bad_transcripts) { next }
        
        # get consequence and severity rank in the current transcript
        consequence = transcript$consequence_terms
        temp_severity = min(severity$rank[severity$consequence %in% consequence])
        consequence = severity$consequence[severity$rank == temp_severity]
        
        # check if this is the most severe consequence; prefer HGNC transcripts
        if (is.na(best_severity) | temp_severity < best_severity |
                (symbol != "HGNC" & transcript$gene_symbol_source == "HGNC" &
                temp_severity == best_severity)) {
            best_severity = temp_severity
            best_transcript = transcript
            symbol = transcript$gene_symbol_source
        }
    }
    
    # sometimes we don't have any useable transcripts, in which case, select
    # any of the transcripts
    if (length(best_transcript) == 1) {
        if (is.na(best_transcript)) {
            best_transcript = find_most_severe_transcript(ensembl_json, exclude_bad=FALSE)
        }
    }
    
    return(best_transcript)
}

#' find the coding length for a ensembl transcript
#'
#' @param transcript_id json data for variant from Ensembl
#' @param build genome build to find consequences on
#' @param verbose flag indicating whether to print variants as they are checked
#'
#' @export
#' @return integer for CDS length.
get_protein_length <-function(transcript_id, build="grch37", verbose=FALSE) {
    ext = "/sequence/id/"
    ext = paste(ext, transcript_id, sep="")
    
    json = try(request_from_ensembl(ext, headers="type=protein", build), silent=TRUE)
    if (class(json) == "try-error") { return(0) }
    json = rjson::fromJSON(json)
    
    return(nchar(json$seq))
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
    
    # make sure the end position is suitable for the Ensembl REST API request
    if (as.numeric(variant[["end_pos"]]) < as.numeric(variant[["start_pos"]])) {
        variant[["end_pos"]] = variant[["start_pos"]]
    }
    
    # define parts of the URL
    ext = "/overlap/region/human/"
    ext = paste(ext, variant[["chrom"]], ":", variant[["start_pos"]], ":", variant[["end_pos"]], "/", sep="")
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]],
            "    ", ext, sep=""))
    }
    
    json = request_from_ensembl(ext, headers="feature=gene", build)
    json = rjson::fromJSON(json)
    
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
    
    # make sure the end position is suitable for the Ensembl REST API request
    if (as.numeric(variant[["end_pos"]]) < as.numeric(variant[["start_pos"]])) {
        variant[["end_pos"]] = variant[["start_pos"]]
    }
    
    # define the URL
    ext = "/sequence/region/human/"
    ext = paste(ext, variant[["chrom"]], ":", variant[["start_pos"]], ":",
        variant[["end_pos"]], ":1", sep="")
    
    if (verbose) {
        print(paste("chr", variant[["chrom"]], ":", variant[["start_pos"]],
            "    ", ext, sep=""))
    }
    
    json = request_from_ensembl(ext, build)
    json = rjson::fromJSON(json)
    
    # return blank string for variants not in genes
    if (length(json) > 0) {
        return(json$seq)
    }
    
    return("")
}
