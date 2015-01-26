# function to open de novo data from the DDD study
#

#' find diagnosed probands in the DDD study, to exclude them from our data
#'
#' @param path path to file defining diagnosed probands
#'
#' @export
#' @return A list containing vectors with DDD IDs, and sex of the diagnosed
#'     probands
get_ddd_diagnosed <- function(path) {
    
    # the diagnoses file is a sheet of an excel file, one row per proband,
    # containing decipher_id, ddd_id, sex, then 11 columns of how many diagnostic
    # variants each proband has been identified with from different mutation
    # categories, eg DNM_SNV (de novo mutation single nucleotide variant),
    # DNM_CNV, AD_INH_SNV (autosomally dominantly inherited SNV).
    
    # read in samples that have been diagnosed, so as to remove from our data
    diagnoses = gdata::read.xls(path, sheet="1133 Diagnoses", stringsAsFactors=FALSE)
    
    # remove samples without DDD IDs
    diagnoses = diagnoses[diagnoses$DDD_ID != "Not in 1133", ]
    
    # find the samples with a diagnosis
    index = rowSums(diagnoses[, c(4:14)]) > 0
    
    diagnosed = list()
    diagnosed$id = diagnoses$DDD_ID[index]
    diagnosed$sex = diagnoses$Sex[index]
    
    return(diagnosed)
}

#' find probands likely to have diagnoses, to exclude them from our data
#'
#' We assume DDD samples that have
#'
#' @param path path to file defining diagnosed probands
#' @export
#'
#' @return A list containing vectors with DDD IDs, and sex of the diagnosed
#'     probands
get_likely_diagnosed <- function(path) {
    
    first_set = get_ddd_diagnosed(path)
    
    DATAFREEZE = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/"
    INDIVIDUALS = file.path(DATAFREEZE, "family_relationships.txt")
    KNOWN_GENES = file.path(DATAFREEZE, "DDG2P_freeze_with_gencode19_genomic_coordinates_20141118_fixed.txt")
    
    probands = read.table(INDIVIDUALS, header=TRUE, stringsAsFactors=FALSE)
    probands = probands[probands$dad_id != 0, ]
    vcf_paths = get_clinical_filtering_vcf_paths(probands$individual_id)
    
    variants = get_variants(vcf_paths)
    
    # find the de novo variants
    de_novos = variants[grepl("deNovo", variants$format_values), ]
    
    # find which of the clinically filtered de novos are also in the set of
    # high quality DDD de novos
    de_novo_key = paste(de_novos$CHROM, de_novos$POS)
    ddd_key = paste(mupit::ddd_de_novos$chrom, mupit::ddd_de_novos$start_pos)
    high_quality = de_novo_key %in% ddd_key
    
    # get the sample IDs and sex for the probands with high quality, clinically
    # filtered de novos
    sample_ids = unique(c(de_novos$sample_id[high_quality], first_set$id))
    sex = probands$sex[probands$individual_id %in% sample_ids]
    
    # If the length of the sex vector is different from the length of the
    # sample_ids vector, then something has ghone wrong, and we will have
    # incorrect stats (particularly when we estimate the number of expected de
    # novos). Stop and figure out the discrepancy before running any more
    # analyses. I know of one discrepant individual in the DDD 4300 trios where
    # the proband had a diagnosis in the 1133 trios, but the proband isn't part
    # of the 4300 trios.
    stopifnot(abs(length(sample_ids) - length(sex)) > 1)
    
    diagnosed = list(id=sample_ids, sex=sex)
    
    return(diagnosed)
}

#' find paths to clinical filtering output for DDD samples
#'
#' @param proband_ids vector of proband IDs
#' @export
#'
#' @return vector of VCF paths
get_clinical_filtering_vcf_paths <- function(proband_ids) {
    DIR = "/lustre/scratch113/projects/ddd/users/ddd/ddd_data_releases/2014-11-04"
    
    # report progress when finding the VCF paths, ~10s for 4300 probands.
    cat("finding DDD VCF paths\n")
    x = 0
    pb = txtProgressBar(min=0, max=length(proband_ids), style=3)
    
    vcf_paths = c()
    for (id in proband_ids) {
        split_path = paste(substring(id, seq(1, nchar(id), 2), seq(2, nchar(id), 2)), collapse="/")
        
        folder = file.path(DIR, split_path, id, "clin_filt")
        
        # often there will be multiple VCFs for a proband, one for each time the
        # clinical filtering has been run for that sample. We simply find the
        # most recent VCF for each proband, by using the dates encoded in the path.
        if (file.exists(folder)) {
            paths = Sys.glob(file.path(folder, paste(id, "clin_filt", "*", "vcf.gz", sep=".")))
            dates = basename(paths)
            dates = sapply(strsplit(dates, "clin_filt."), "[", 2)
            dates = sapply(strsplit(dates, ".vcf.gz"), "[", 1)
            
            vcf_path = paths[order(as.Date(dates, format="%Y-%m-%d"), decreasing=TRUE)[1]]
            vcf_paths = c(vcf_paths, vcf_path)
        }
        x = x + 1
        if (x %% 20 == 0) { setTxtProgressBar(pb, x) }
    }
    
    close(pb)
    
    return(vcf_paths)
}

#' get standardised de novo data for DDD study.
#'
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
standardise_ddd_de_novos <- function() {
    
    # load a set of DDD de novos, that are the independent event (ie duplicate
    # de novos within families have been removed).
    variants = read.table(file.path("data-raw", "de_novo_datasets",
        "de_novos.ddd_4k.ddd_only.txt"), header=TRUE,
        sep="\t", stringsAsFactors=FALSE, comment.char="")
    
    # standardise the SNV or INDEL flag
    variants$type = "indel"
    variants$type[variants$var_type == "DENOVO-SNP"] = "snv"
    
    # standardise the columns, and column names
    variants$person_id = gsub(" ", "", variants$person_stable_id)
    variants$chrom = gsub(" ", "", variants$chrom)
    variants$start_pos = gsub(" ", "", variants$pos)
    variants$ref_allele = gsub(" ", "", variants$ref)
    variants$alt_allele = gsub(" ", "", variants$alt)
    variants$end_pos = as.character(as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1)
    
    variants$hgnc = variants$symbol
    # vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    # variants$consequence = sapply(vep, "[", 1)
    # variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "ddd_unpublished"
    variants$publication_doi = NA
    variants$study_phenotype = "developmental_disorders"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype", "type"))
    
    return(variants)
}

#' get de novo data for DDD study.
#'
#' @param diagnosed dataframe of samples in DDD with a diagnosis, who should be
#'     excluded from analyses
#' @param subset vector of sample IDs, if we want to restrict to a specific
#'         set of probands, such as when we are analysing a phenotype-specific
#'         subset of the DDD.
#'
#' @export
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
get_ddd_de_novos <- function(diagnosed=NULL, subset=NULL) {
    
    # try and use the data from the package, if it is available, otherwise
    # generate the dataset
    if (exists("ddd_de_novos") & nrow(ddd_de_novos) > 0) {
        variants = ddd_de_novos
     } else {
        variants = standardise_ddd_de_novos()
        ddd_de_novos = variants
        save(ddd_de_novos, file="data/ddd_de_novos.rda", compress="xz")
    }
    
    # remove diagnosed patients, if maximising power
    if (!is.null(diagnosed)[1]) {
        variants = variants[!(variants$person_id %in% diagnosed$id), ]
    }
    
    # sometimes we only want to use a subset of the DDD, such as when we are
    # investigating a single disease type, e.g. seizures. Then we will have
    # specified the DDD sample IDs that we wish to restrict to.
    if (!is.null(subset)[1]) {
        variants = variants[variants$person_id %in% subset, ]
    }
    
    return(variants)
}
