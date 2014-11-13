# Script to generate processed datasets for the mupit package

library(mupit)

#' get de novo data for Rauch et al. intellectual disability exome study
#' 
#' De novo mutation data sourced from supplementary tables 2 and 3 from
#' Rauch et al. (2012) Lancet 380:1674-1682 
#' doi: 10.1016/S0140-6736(12)61480-9
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
prepare_rauch_de_novos <- function() {
   
    rauch_de_novos = read.delim(file.path("data-raw", "de_novo_datasets", "rauch_v2.txt"), header=TRUE, colClasses = "character")
    
    # standardise the columns, and column names
    rauch_de_novos = subset(rauch_de_novos, select = c("INFO.HGNC", "INFO.CQ", "POS", "CHROM", "TYPE"))
    names(rauch_de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    rauch_de_novos$study = "rauch"
    
    save(rauch_de_novos, file="data/rauch_de_novos.rda")
}

#' get de novo data for De Ligt et al. intellectual disability exome study
#' 
#' De novo mutation data sourced from supplementary table 3 from
#' De Ligt et al. (2012) N Engl J Med (2012) 367:1921-1929
#' doi: 10.1056/NEJMoa1206524
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
prepare_deligt_de_novos <- function() {
    
    deligt_de_novos = read.delim(file.path("data-raw", "de_novo_datasets", "deligt_v2.txt"), header=TRUE, colClasses = "character")
    
    # standardise the columns, and column names
    deligt_de_novos = subset(deligt_de_novos, select = c("INFO.HGNC", "INFO.CQ", "POS", "CHROM", "TYPE"))
    names(deligt_de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    deligt_de_novos$study = "deligt"
    
    save(deligt_de_novos, file="data/deligt_de_novos.rda")
}

#' get de novo data for the Epi4K epilepsy exome study
#' 
#' De novo mutation data sourced from supplementary table 2 from
#' Allen et al. (2013) Nature 501:217-221 
#' doi: 10.1038/nature12439
#' 
#' TODO: swap out this data for that used in the more recent EPI4K publication:
#' TODO: American Journal of Human Genetics (2014) 95:360-370
#' TODO: doi: 10.1016/j.ajhg.2014.08.013
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
prepare_epi4k_de_novos <- function() {
    
    epi4k_de_novos = gdata::read.xls("http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0002929714003838/1-s2.0-S0002929714003838-mmc2.xlsx/276895/FULL/S0002929714003838/bf21945d72e3297fc44969dc0296f4f1/mmc2.xlsx", stringsAsFactors=FALSE)
    
    # the excel table contains three final tail rows which do not contain
    # tabular data, so we remove these
    epi4k_de_novos = epi4k_de_novos[1:(nrow(epi4k_de_novos) - 3), ]
    
    epi4k_de_novos = prepare_coordinates_with_allele(epi4k_de_novos, 
        "hg19.coordinates..chr.position.", "Ref.Alt.alleles")
    
    # get the HGNC symbol
    # NOTE: the variant with the most severe consequence might not necessarily  
    # NOTE: be within the listed gene. I haven't accounted for this yet.
    epi4k_de_novos$consequence = apply(epi4k_de_novos, 1, get_most_severe_vep_consequence, verbose=TRUE)
    
    # get the hgnc symbol, and clean any anomalies
    epi4k_de_novos$hgnc = epi4k_de_novos$Gene
    epi4k_de_novos$hgnc = gsub(" \\(MLL4\\)", "", epi4k_de_novos$hgnc)
    epi4k_de_novos$hgnc = gsub(" \\^", "", epi4k_de_novos$hgnc)
    
    # set the SNV/INDEL type
    epi4k_de_novos$type = "INDEL"
    epi4k_de_novos$type[epi4k_de_novos$start_pos == epi4k_de_novos$end_pos] = "SNV"
    
    epi4k_de_novos = subset(epi4k_de_novos, select = c("hgnc", "consequence", "start_pos", "chrom", "type"))
    names(epi4k_de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    epi4k_de_novos$study = "epi4k"
    
    save(epi4k_de_novos, file="data/epi4k_de_novos.rda")
}

#' get de novo data from the Sanders et al autism exome study
#' 
#' Supplementary table 2 (where the excel sheets for the probands and 
#' siblings have been combined) from:
#' Sanders et al. (2012) Nature 485:237-241 
#' doi: 10.1038/nature10945
#' 
#' @return data frame of de novos, with standardised genome coordinates and VEP
#      consequences for each variant
prepare_sanders_de_novos <- function() {
    url = "http://www.nature.com/nature/journal/v485/n7397/extref/nature10945-s3.xls"
    sanders_probands = gdata::read.xls(url, sheet="Probands", stringsAsFactors=FALSE)
    sanders_siblings = gdata::read.xls(url, sheet="Siblings", stringsAsFactors=FALSE)
    variants = rbind(sanders_probands, sanders_siblings)
    
    variants$chrom = gsub("chr", "", variants$Chr.1)
    variants$start_pos = gsub(" ", "", variants$Pos..hg19.)
    variants$ref_allele = variants$Ref
    variants$alt_allele = variants$Alt
    
    variants$temp = NA
    
    indels = grep(":", variants$alt_allele)
    
    variants$ref_allele[indels] = sapply(strsplit(variants$alt_allele[indels], ":"), "[[", 2)
    variants$alt_allele[indels] = "-"
    
    variants$end_pos = variants$start_pos
    variants$end_pos[indels] = as.numeric(variants$end_pos[indels]) + nchar(variants$ref_allele[indels])
    
    variants = fix_het_alleles(variants)
    
    variants$consequence = NA
    for (row_num in 1:nrow(variants)) {
        variants$consequence[row_num] = get_most_severe_vep_consequence(variants[row_num, ], verbose=TRUE)
    }
    variants$hgnc = variants$Gene
    
    return(variants)
}

#' get de novo data from the O'Roak et al autism exome study
#' 
#' Supplementary table 3 from:
#' O'Roak et al. (2012) Nature 485:246-250
#' doi: 10.1038/nature10989
#' 
#' @return data frame of de novos, with standardised genome coordinates and VEP
#      consequences for each variant
prepare_oroak_de_novos <- function() {
    url = "http://www.nature.com/nature/journal/v485/n7397/extref/nature10989-s2.xls"
    variants = gdata::read.xls(url, sheet="Supplementary Table 3", stringsAsFactors=FALSE)
    
    # the excel table contains two final tail rows which do not contain
    # tabular data, so we remove these
    variants = variants[1:(nrow(variants) - 2), ]
    
    # standardise the chrom, position and allele column names
    variants$chrom = variants$Chromosome
    variants$start_pos = variants$Position..hg19.
    variants$ref_allele = variants$Ref
    variants$alt_allele = variants$Allele
    
    # sort out the "complex" allele events
    variants$alt_allele[61] = "S"
    variants$alt_allele[78] = "1D, -G"
    variants$alt_allele[116] = "1I, +C"
    variants$alt_allele[133] = "R"
    variants$alt_allele[185] = "1D, -G"
    
    # fix the alleles and positions for insertions and deletions
    deletions = grep("D, -", variants$alt_allele)
    insertions = grep("I, +", variants$alt_allele)
    
    variants$ref_allele[deletions] = sapply(strsplit(variants$alt_allele[deletions], "D, -"), "[[", 2)
    variants$alt_allele[deletions] = "-"
    variants$ref_allele[insertions] = "N"
    variants$alt_allele[insertions] = paste("A", sapply(strsplit(variants$alt_allele[insertions], "I, \\+"), "[[", 2), sep="")
    
    # get the end coordinate, including those for the insertions and deletions
    variants$end_pos = variants$start_pos
    variants$end_pos[deletions] = as.numeric(variants$end_pos[deletions]) + nchar(variants$ref_allele[deletions]) - 1
    
    variants = fix_het_alleles(variants)
    variants$consequence = apply(variants, 1, get_most_severe_vep_consequence, verbose=TRUE)
    variants$hgnc = variants$Gene_SeaSeq
    
    return(variants)
}

#' get de novo data from the 2012 Iossifov et al autism exome study in Neuron
#' 
#' Supplementary table 1 (where the non-coding SNVs have been excluded) from:
#' Iossifov et al. (2012) Neuron 74:285-299
#' doi: 10.1016/j.neuron.2012.04.009
#' 
#' @return data frame of de novos, with standardised genome coordinates and VEP
#      consequences for each variant
prepare_iossifov_neuron_de_novos <- function() {
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0896627312003406/1-s2.0-S0896627312003406-mmc2.xlsx/272195/FULL/S0896627312003406/26c5ba3b72a2410ef43fec52a40f35e6/mmc2.xlsx"
    variants = gdata::read.xls(url, sheet="SNV.v4.1-normlized", stringsAsFactors=FALSE)
    
    # trim out the low quality de novos (as defined by a flag in the table)
    variants = variants[variants$SNVFilter == 1, ]
    
    # get the coordinates and VEP consequence
    variants = prepare_coordinates_with_allele(variants, "location", "variant")
    variants$consequence = apply(variants, 1, get_most_severe_vep_consequence, verbose=TRUE)
    
    # get hgnc symbols for the few genes that have missing values
    no_gene = variants$effectGenes == ""
    variants$effectGenes[no_gene] = apply(variants[no_gene, ], 1, get_gene_id_for_variant)
    
    # extract the hgnc symbol from a gene:consequence string
    variants$hgnc = sapply(strsplit(variants$effectGenes, ":"), "[[", 1)
    
    # exclude de novos not located within the coding sequence of a gene
    noncoding_consequences = c("non_coding_transcript_exon_variant",
        "3_prime_UTR_variant", "intron_variant", "downstream_gene_variant",
        "5_prime_UTR_variant", "upstream_gene_variant")
    variants = variants[!(variants$consequence %in% noncoding_consequences), ]
    
    return(variants)
}

#' get de novo data from the 2014 Iossifov et al. autism exome study in Nature
#' 
#' De novo mutation data sourced from Supplementary table 2:
#' Iossifov et al. (2013) Nature 498:220-223
#' doi: 10.1038/nature13908
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
prepare_iossifov_nature_de_novos <-function() {
    tmpdir = tempdir()
    path = tempfile(tmpdir=tmpdir)
    
    # obtain the dataframe of de novo variants
    download.file("http://www.nature.com/nature/journal/vaop/ncurrent/extref/nature13908-s2.zip", path)
    unzip(path, files="nature13908-s2/Supplementary Table 2.xlsx", exdir=tmpdir)
    variants = gdata::read.xls(file.path(tmpdir, "nature13908-s2", "Supplementary Table 2.xlsx"), stringsAsFactors=FALSE)
    
    variants = prepare_coordinates(variants, "location", "vcfVariant")
    
    # NOTE: the variant with the most severe consequence might not necessarily  
    # NOTE: be within the listed gene. I haven't accounted for this yet.
    variants$consequence = apply(variants, 1, get_most_severe_vep_consequence, verbose=TRUE)
    
    # get the HGNC symbol
    no_gene = variants$effectGene == ""
    variants$effectGene[no_gene] = apply(variants[no_gene, ], 1, get_gene_id_for_variant)
    variants$hgnc = variants$effectGene
    
    # set the SNV/INDEL type
    variants$type = "INDEL"
    variants$type[variants$start_pos == variants$end_pos] = "SNV"
    
    unlink(path)
    
    # remove the de novos identified in the Iossifov publication, since they are
    # also present in the autism_de_novos dataset in this package.
    # variants = variants[variants$IossifovWE2012 != "yes", ]
    
    return(variants)
    
    # variants = subset(variants, select = c("hgnc", "consequence", "start_pos", "chrom", "type"))
    # names(variants) = c("hgnc", "consequence", "position", "chrom", "type")
    # variants$study = "iossifov"
    
    # save(iossifov_de_novos, file="data/iossifov_de_novos.rda")
}

#' get de novo data from autism exome studies
#' 
#' I think this data was obtained from studies published using the Simon's
#' Simplex Collection data (http://sfari.org/resources/simons-simplex-collection)
#' 
#' Supplementary table 2 (where the excel sheets for the probands and 
#' siblings have been combined) from:
#' Sanders et al. (2012) Nature 485:237-241 
#' doi: 10.1038/nature10945
#' 
#' Supplementary table 3 from:
#' O'Roak et al. (2012) Nature 485:246-250
#' doi: 10.1038/nature10989
#' 
#' Supplementary table 1 (where the non-coding SNVs have been excluded) from:
#' Iossifov et al. (2012) Neuron 74:285-299
#' doi: 10.1016/j.neuron.2012.04.009
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
prepare_autism_de_novos <- function() {
    
    sanders = prepare_sanders_de_novos()
    oroak = prepare_oroak_de_novos()
    iossifov_neuron = prepare_iossifov_neuron_de_novos()
    iossifov_nature = prepare_iossifov_nature_de_novos()
    
    autism_de_novos = read.delim(file.path("data-raw", "de_novo_datasets", "autism_v3_PJ.txt"), header=TRUE, colClasses = "character")
    
    # select only de novos in probands
    autism_de_novos = autism_de_novos[which(autism_de_novos$pheno == "Pro"), ]
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(nchar(autism_de_novos$ref.1) != nchar(autism_de_novos$var))
    autism_de_novos$TYPE = "INDEL"
    autism_de_novos$TYPE[TYPE.index] = "SNV"
    
    # standardise the columns, and column names
    autism_de_novos = subset(autism_de_novos, select = c("INFO.HGNC", "INFO.CQ", "pos", "CHROM", "TYPE"))
    names(autism_de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    autism_de_novos$study = "autism"
    
    save(autism_de_novos, file="data/autism_de_novos.rda")
}

#' get de novo data from Fromer et al. schizophrenia exome study
#' 
#' De novo mutation data sourced from Supplementary table 1:
#' Fromer et al. (2014) Nature 506:179-184
#' doi: 10.1038/nature12929
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
prepare_fromer_de_novos <- function() {
    
    fromer_de_novos = read.delim(file.path("data-raw", "de_novo_datasets", "fromer_v2.txt"), header=TRUE, colClasses = "character")
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(abs(nchar(fromer_de_novos$Reference.allele) - nchar(fromer_de_novos$Alternate.allele)) == 0)
    fromer_de_novos$TYPE = "INDEL"
    fromer_de_novos$TYPE[TYPE.index] = "SNV"
    
    # standardise the columns, and column names
    fromer_de_novos = subset(fromer_de_novos, select = c("INFO.HGNC", "INFO.CQ", "pos", "chrom", "TYPE"))
    names(fromer_de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    fromer_de_novos$study = "fromer"
    
    save(fromer_de_novos, file="data/fromer_de_novos.rda")
}

#' get de novo data from Zaidi et al. congenital heart disease exome study
#' 
#' De novo mutation data sourced from Supplementary table 4:
#' Zaidi et al. (2013) Nature 498:220-223
#' doi: 10.1038/nature12141
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
prepare_zaidi_de_novos <- function() {
    
    # could only include syndromic DNMs
    zaidi_de_novos = read.delim(file.path("data-raw", "de_novo_datasets", "zaidi_VEP.txt"), header=TRUE, colClasses = "character")
    
    # remove DNMs in controls
    zaidi_de_novos = zaidi_de_novos[-which(zaidi_de_novos$Primary_Cardiac_Class == "Control"), ]
    
    # standardise the SNV or INDEL flag
    TYPE.index = which(abs(nchar(zaidi_de_novos$ref) - nchar(zaidi_de_novos$alt)) != 0)
    TYPE.index = sort(unique(c(TYPE.index, which(zaidi_de_novos$ref == "-"), which(zaidi_de_novos$alt == "-"))))
    zaidi_de_novos$TYPE = "SNV"
    zaidi_de_novos$TYPE[TYPE.index] = "INDEL"
    
    # standardise the columns, and column names
    zaidi_de_novos = subset(zaidi_de_novos, select = c("INFO.HGNC", "INFO.CQ", "pos", "chrom", "TYPE"))
    names(zaidi_de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    zaidi_de_novos$study = "zaidi"
    
    save(zaidi_de_novos, file="data/zaidi_de_novos.rda")
}

prepare_rauch_de_novos()
prepare_deligt_de_novos()
prepare_epi4k_de_novos()
prepare_autism_de_novos()
prepare_fromer_de_novos()
prepare_zaidi_de_novos()
prepare_iossifov_de_novos()
