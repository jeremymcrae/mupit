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

#' get de novo data from Iossifov et al. congenital heart disease exome study
#' 
#' De novo mutation data sourced from Supplementary table 2:
#' Iossifov et al. (2013) Nature 498:220-223
#' doi: 10.1038/nature13908
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
prepare_iossifov_de_novos <-function() {
    tmpdir = tempdir()
    path = tempfile(tmpdir=tmpdir)
    
    # obtain the dataframe of de novo variants
    download.file("http://www.nature.com/nature/journal/vaop/ncurrent/extref/nature13908-s2.zip", path)
    unzip(path, files="nature13908-s2/Supplementary Table 2.xlsx", exdir=tmpdir)
    iossifov_de_novos = gdata::read.xls(file.path(tmpdir, "nature13908-s2", "Supplementary Table 2.xlsx"), stringsAsFactors=FALSE)
    
    iossifov_de_novos = prepare_coordinates(iossifov_de_novos, "location", "vcfVariant")
    
    # get the HGNC symbol
    # NOTE: the variant with the most severe consequence might not necessarily  
    # NOTE: be within the listed gene. I haven't accounted for this yet.
    iossifov_de_novos$hgnc = iossifov_de_novos$effectGene
    iossifov_de_novos$consequence = apply(iossifov_de_novos, 1, get_most_severe_vep_consequence, verbose=TRUE)
    
    # set the SNV/INDEL type
    iossifov_de_novos$type = "INDEL"
    iossifov_de_novos$type[iossifov_de_novos$start_pos == iossifov_de_novos$end_pos] = "SNV"
    
    unlink(path)
    
    # remove the de novos identified in the Iossifov publication, since they are
    # also present in the autism_de_novos dataset in this package.
    iossifov_de_novos = iossifov_de_novos[iossifov_de_novos$IossifovWE2012 != "yes", ]
    
    iossifov_de_novos = subset(iossifov_de_novos, select = c("hgnc", "consequence", "start_pos", "chrom", "type"))
    names(iossifov_de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    iossifov_de_novos$study = "iossifov"
    
    save(iossifov_de_novos, file="data/iossifov_de_novos.rda")
}

prepare_rauch_de_novos()
prepare_deligt_de_novos()
prepare_epi4k_de_novos()
prepare_autism_de_novos()
prepare_fromer_de_novos()
prepare_zaidi_de_novos()
prepare_iossifov_de_novos()
