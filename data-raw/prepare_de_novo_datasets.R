# Script to generate processed datasets for the mupit package

library(mupit)

# exclude de novos not located within the coding sequence of a gene
NONCODING_CONSEQUENCES = c("non_coding_transcript_exon_variant",
    "3_prime_UTR_variant", "intron_variant", "downstream_gene_variant",
    "5_prime_UTR_variant", "upstream_gene_variant", "intergenic_variant",
    "regulatory_region_variant")

#' get de novo data for Rauch et al. intellectual disability exome study
#'
#' De novo mutation data sourced from supplementary tables 2 and 3 from
#' Rauch et al. (2012) Lancet 380:1674-1682
#' doi: 10.1016/S0140-6736(12)61480-9
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
rauch_de_novos <- function() {
    
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0140673612614809/1-s2.0-S0140673612614809-mmc1.pdf/271074/FULL/S0140673612614809/55b26043f4a279334b3a5ec00b9faf4b/mmc1.pdf"
    
    # obtain the supplementary material
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    # get the lines corresponding to each table
    table_s1_text = test$content[333:529]
    table_s2_text = test$content[769:847]
    table_s3_text = test$content[879:890]
    
    # clean up table S1
    table_s1_text = gsub("^[ \t\f]+", "", table_s1_text) # trim leading whitespace
    split_strings = strsplit(table_s1_text, "[ \t]{2,}")
    split_strings = split_strings[sapply(split_strings, length) > 1]
    split_strings = lapply(split_strings, function(x) x[1:3])
    table_s1 = data.frame(t(data.frame(split_strings)))
    row.names(table_s1) = 1:nrow(table_s1)
    
    # only select the lines with a single-letter sex code
    use_lines = apply(table_s1, 1, function(x) x[2] %in% c("f", "m") | x[3] %in% c("f", "m"))
    table_s1 = table_s1[use_lines, ]
    table_s1$sex = apply(table_s1, 1, function(x) x[x == "f" | x == "m"])
    names(table_s1) = c("person_id", "column2", "column3", "sex")
    table_s1 = table_s1[, c("person_id", "sex")]
    table_s1$sex = gsub("m", "male", table_s1$sex)
    table_s1$sex = gsub("f", "female", table_s1$sex)
    table_s1$person_id = gsub("-", "", table_s1$person_id)
    table_s1$person_id = as.character(table_s1$person_id)
    table_s1$person_id[table_s1$person_id == "TUTLN112014"] = "TUTLN"
    
    # clean up table S2
    table_s2_text = gsub("^[ \t\f]+", "", table_s2_text) # trim leading whitespace
    split_strings = strsplit(table_s2_text, "[ \t]+")
    split_strings = iconv(unlist(split_strings), "latin1", "ASCII", sub="")
    table_s2 = as.data.frame(matrix(split_strings, ncol=12, byrow=TRUE))
    names(table_s2) = c("person_id", "hgnc", "type", "hgvs_genomic")
    
    # clean up table S3
    table_s3_text = gsub("^[ \t\f]+", "", table_s3_text) # trim leading whitespace
    split_strings = strsplit(table_s3_text, "[ \t]+")
    split_strings[7][[1]] = c(split_strings[[7]], "") # one row lacks an entry
    split_strings = iconv(unlist(split_strings), "latin1", "ASCII", sub="")
    table_s3 = as.data.frame(matrix(split_strings, ncol=9, byrow=TRUE))
    names(table_s3) = c("person_id", "hgnc", "hgvs_genomic")
    table_s3$type = "synonymous"
    
    # standardise and merge the two tables
    table_s2 = subset(table_s2, select=c(person_id, hgnc, type, hgvs_genomic))
    table_s3 = subset(table_s3, select=c(person_id, hgnc, type, hgvs_genomic))
    variants = rbind(table_s2, table_s3)
    
    variants = fix_coordinates_with_hgvs_genomic(variants, "hgvs_genomic")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    # include the sex codes for each de novo
    variants = merge(variants, table_s1, by="person_id", all.x=TRUE)
    
    # define the study details
    variants$study_code = "rauch_lancet_2012"
    variants$publication_doi = "10.1016/S0140-6736(12)61480-9"
    variants$study_phenotype = "intellectual_disability"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}

#' get de novo data for De Ligt et al. intellectual disability exome study
#'
#' De novo mutation data sourced from supplementary table 3 from
#' De Ligt et al. (2012) N Engl J Med 367:1921-1929
#' doi: 10.1056/NEJMoa1206524
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
deligt_de_novos <- function() {
    
    url = "http://www.nejm.org/doi/suppl/10.1056/NEJMoa1206524/suppl_file/nejmoa126524_appendix.pdf"
    
    temp = "temp.html"
    cookie = "cookie.txt"
    system(paste("wget --cookies=on --keep-session-cookies --save-cookies=", cookie, " ", url, " -O", temp, sep=""))
    
    # obtain the supplementary material
    path = "corrupted_pdf.pdf"
    system(paste("wget --referer=", url, " --cookies=on --load-cookies=", cookie, " --keep-session-cookies --save-cookies=", cookie, " ", url, " -O ", path, sep=""))
    
    # repair the pdf with ghostscript
    repaired = "repaired.pdf"
    system(paste("gs -o ", repaired, " -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress ", path, sep=""))
    
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri=repaired), language = "en", id = "id1")
    
    # delete the pdf and temporary files
    unlink(temp)
    unlink(cookie)
    unlink(path)
    unlink(repaired)
    
    # find the sex of each proband
    text = test$content[174:1160]
    text = gsub("^[ \t\f]+", "", text)
    gender = data.frame(text=text[grepl("^Trio", gender)], sex=NA)
    gender$text = as.character(gender$text)
    gender$sex[grepl("[Tt]his (boy|male)+", gender$text)] = "male"
    gender$sex[grepl("[Tt]his (girl|female)+", gender$text)] = "female"
    gender$person_id = sapply(strsplit(gender$text, "( |-)+"), "[", 2)
    gender$sex[gender$person_id == "69"] = "female"
    
    # get the lines corresponding to each table
    table_s3_text = test$content[1709:1795]
    
    # clean up table S2
    table_s3_text = gsub("^[ \t\f]+", "", table_s3_text) # trim leading whitespace
    table_s3_text = table_s3_text[table_s3_text != ""]
    table_s3_text = table_s3_text[table_s3_text != "#"]
    table_s3_text = table_s3_text[lapply(table_s3_text, nchar) > 5]
    split_strings = strsplit(table_s3_text, "[ \t]+")
    
    # cull it down to the first few entries in each line
    variants = data.frame(t(sapply(split_strings, "[", 1:6)))
    names(variants) = c("person_id", "hgnc", "hgvs_genomic", "transcript", "hgvs_cdna", "hgvs_protein")
    
    # fix the hgvs genomic string
    variants$hgvs_genomic = gsub("\\(GRCh37\\)g", ":g", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("\\(GRCh37\\)", "", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("\\(GRCH37\\)", "", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("Chr", "chr", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("-", "_", variants$hgvs_genomic)
    
    # convert a position with NCBI36 genome assembly coordinates
    variants$hgvs_genomic[27] = "chr19:g.53958839G>C"
    
    variants = fix_coordinates_with_hgvs_genomic(variants, "hgvs_genomic")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants = merge(variants, gender, by="person_id")
    
    variants$study_code = "deligt_nejm_2012"
    variants$publication_doi = "10.1056/NEJMoa1206524"
    variants$study_phenotype = "intellectual_disability"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(list(variants=variants, gender=gender))
}

#' get de novo data for Gilissen et al. intellectual disability exome study
#'
#' De novo mutation data sourced from supplementary table 8 from:
#' Gilissen et al. (2014) Nature 511:344-347
#' doi: 10.1038/nature13394
#'
#' NOTE: the number of males and females will form part of the De Ligt counts,
#' NOTE: since the Gilissen samples were in the De Ligt study, they just didn't
#' receive a diagnosis in the original study
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
gilissen_de_novos <- function(deligt, gender) {
    
    url = "http://www.nature.com/nature/journal/v511/n7509/extref/nature13394-s1.pdf"
    
    # obtain the supplementary material
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    # get the lines corresponding to each table
    table_s8 = test$content[1434:1563]
    
    # clean up table S2
    table_s8 = gsub("^[ \t]+", "", table_s8) # trim leading whitespace
    
    # fix the lines that come after the page breaks
    table_s8 = table_s8[table_s8 != ""]
    table_s8 = table_s8[table_s8 != "#"]
    table_s8 = table_s8[!grepl("NATURE.COM", table_s8)]
    table_s8 = table_s8[!grepl("SUPPLEMENTARY INFORMATION", table_s8)]
    
    # one line in the table is spread across three lines, rather than fixing
    # this in code, just re-insert the correct line
    table_s8[50] = "25  SATB2     Chr2(GRCh37):g.200213667_200213668insGTTGCCTTACAA    NM_001172517.1:c.929_930insTTGTAAGGCAAC  p.Gln310delinsHis_CysLys_AlaThr  Insertion  no  K  F  D  C  B  G  M  Known"
    
    # trim the too short lines (which only contain page numbers)
    table_s8 = table_s8[lapply(table_s8, nchar) > 5]
    
    # split the table, and drop the erroneous lines from the bad line
    split_strings = strsplit(table_s8, "[ \t]+")
    split_strings[50] = NULL
    split_strings[49] = NULL
    
    # cull it down to the first few entries in each line
    variants = data.frame(t(sapply(split_strings, "[", 1:6)))
    names(variants) = c("person_id", "hgnc", "hgvs_genomic", "hgvs_transcript", "hgvs_protein", "type")
    
    # fix the hgvs genomic string
    variants$hgvs_genomic = gsub("\\(GRCh37\\)g", ":g", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("\\(GRCh37\\)", "", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("\\(GRCH37\\)", "", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("Chr", "chr", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("-", "_", variants$hgvs_genomic)
    
    variants = fix_coordinates_with_hgvs_genomic(variants, "hgvs_genomic")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "gilissen_nature_2014"
    variants$publication_doi = "10.1038/nature13394"
    variants$study_phenotype = "intellectual_disability"
    
    variants = merge(variants, gender, by="person_id", all.x=TRUE)
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    # remove the variants that were present in the De Ligt study
    temp = rbind(deligt, variants)
    temp$dups = duplicated(temp[, 2:7])
    variants = variants[!(temp$dups[temp$study_code == "gilissen_nature_2014"]), ]
    
    return(variants)
}

#' get de novo data for the Epi4K epilepsy exome study
#'
#' De novo mutation data from the most recent EPI4K publication:
#' Supplementary table 1:
#' American Journal of Human Genetics (2014) 95:360-370
#' doi: 10.1016/j.ajhg.2014.08.013
#'
#' This incorporates the de novo mutation data from supplementary table 2 of:
#' Allen et al. (2013) Nature 501:217-221
#' doi: 10.1038/nature12439
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
epi4k_de_novos <- function() {
    
    system("wget https://catalog.coriell.org/0/Excel/NINDS/Epilepsy_affected_clinicaldata.xlsx")
    clinical_data = gdata::read.xls("Epilepsy_affected_clinicaldata.xlsx", sheet="Epilepsy Clinical Data", stringsAsFactors=FALSE)
    clinical_data$sex = NA
    clinical_data$sex[clinical_data$Gender == "F"] = "female"
    clinical_data$sex[clinical_data$Gender == "M"] = "male"
    clinical_data = clinical_data[, c("Catalog_ID", "sex")]
    unlink("Epilepsy_affected_clinicaldata.xlsx")
    
    # obtain the supplementary material
    url = "http://www.nature.com/nature/journal/v501/n7466/extref/nature12439-s1.pdf"
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    table_s1 = test$content[11:175]
    table_s1 = gsub("^[ \t\f]+", "", table_s1)
    table_s1 = table_s1[!grepl("N A T U R E", table_s1)]
    table_s1 = table_s1[!grepl("SUPPLEMENTARY INFORMATION", table_s1)]
    table_s1 = table_s1[!grepl("Trio", table_s1)]
    table_s1 = table_s1[!grepl("Coriell ID", table_s1)]
    table_s1 = table_s1[!grepl("Proband", table_s1)]
    table_s1 = table_s1[nchar(table_s1) > 1]
    table_s1 = strsplit(table_s1, "[ ]+")
    table_s1 = data.frame(t(data.frame(table_s1)))
    row.names(table_s1) = 1:nrow(table_s1)
    part1 = table_s1[, 1:4]
    names(part1) = c("person_id", "proband", "father", "mother")
    part2 = table_s1[, 5:8]
    names(part2) = c("person_id", "proband", "father", "mother")
    
    samples = rbind(part1, part2)
    samples = merge(samples, clinical_data, by.x="proband", by.y="Catalog_ID")
    gender = samples[, c("person_id", "sex")]
    
    # now get the supplementary material from the AJHG paper, for the probands
    # who do not have sex info in the Nature paper supplementary material
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0002929714003838/1-s2.0-S0002929714003838-mmc1.pdf/276895/FULL/S0002929714003838/4a001cfb9c56ad9275203cfbe75e9e9e/mmc1.pdf"
    
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    # clean up the text that forms part of supplementary table S2.
    table_s2 = test$content[200:320]
    table_s2 = gsub("^[ \t\f]+", "", table_s2)
    table_s2 = table_s2[!grepl("^Trio", table_s2)]
    table_s2 = strsplit(table_s2, "[ ]+")
    table_s2 = table_s2[sapply(table_s2, length) > 14]
    table_s2 = sapply(table_s2, function(x) x[c(1, 3)])
    
    table_s2 = data.frame(t(table_s2))
    names(table_s2) = c("person_id", "sex")
    table_s2$person_id = gsub("(LGS|IS)", "", table_s2$person_id)
    table_s2$sex = sapply(strsplit(as.character(table_s2$sex), "/"), "[", 1)
    table_s2$sex[table_s2$sex == "F"] = "female"
    table_s2$sex[table_s2$sex == "M"] = "male"
    
    gender = rbind(gender, table_s2)
    
    variants = gdata::read.xls("http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0002929714003838/1-s2.0-S0002929714003838-mmc2.xlsx/276895/FULL/S0002929714003838/bf21945d72e3297fc44969dc0296f4f1/mmc2.xlsx", stringsAsFactors=FALSE)
    
    # the excel table contains three final tail rows which do not contain
    # tabular data, so we remove these
    variants = variants[1:(nrow(variants) - 3), ]
    
    variants = fix_coordinates_with_allele(variants,
        "hg19.coordinates..chr.position.", "Ref.Alt.alleles")
    
    # get the hgnc symbol, and clean any anomalies
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Child.ID
    variants$study_code = "epi4k_ajhg_2014"
    variants$publication_doi = "10.1016/j.ajhg.2014.08.013"
    variants$study_phenotype = "epilepsy"
    
    variants = merge(variants, gender, by="person_id", all.x=TRUE)
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
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
sanders_de_novos <- function() {
    
    url = "http://www.nature.com/nature/journal/v485/n7397/extref/nature10945-s2.xls"
    samples = gdata::read.xls(url, stringsAsFactors=FALSE)
    samples = samples[!(samples$Role %in% c("Mother", "Father")), ]
    samples$sex = tolower(samples$Gender)
    samples$person_id = samples$Sample
    gender = samples[, c("person_id", "sex")]
    
    url = "http://www.nature.com/nature/journal/v485/n7397/extref/nature10945-s3.xls"
    sanders_probands = gdata::read.xls(url, sheet="Probands", stringsAsFactors=FALSE)
    sanders_siblings = gdata::read.xls(url, sheet="Siblings", stringsAsFactors=FALSE)
    variants = rbind(sanders_probands, sanders_siblings)
    
    variants$chrom = gsub("chr", "", variants$Chr.1)
    variants$start_pos = gsub(" ", "", variants$Pos..hg19.)
    variants$end_pos = variants$start_pos
    variants$ref_allele = variants$Ref
    variants$alt_allele = variants$Alt
    
    # get the correct ref and alt alleles for indels
    indels = grep(":", variants$alt_allele)
    temp_distance = nchar(sapply(strsplit(variants$alt_allele[indels], ":"), "[[", 2))
    variants$alt_allele[indels] = apply(variants[indels, ], 1, get_sequence_in_region)
    variants$end_pos[indels] = as.numeric(variants$end_pos[indels]) + temp_distance
    variants$ref_allele[indels] = apply(variants[indels, ], 1, get_sequence_in_region)
    
    variants = fix_het_alleles(variants)
    
    # get the HGNC symbol and VEP consequence for each variant
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Child_ID
    variants$study_code = "sanders_nature_2012"
    variants$publication_doi = "10.1038/nature10945"
    variants$study_phenotype = "autism"
    
    variants = merge(variants, gender, by="person_id", all.x=TRUE)
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
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
oroak_de_novos <- function() {
    url = "http://www.nature.com/nature/journal/v485/n7397/extref/nature10989-s2.xls"
    variants = gdata::read.xls(url, sheet="Supplementary Table 3", stringsAsFactors=FALSE)
    samples = gdata::read.xls(url, sheet="Supplementary Table 1", stringsAsFactors=FALSE)
    samples$person_id = samples$child
    gender = samples[, c("person_id", "sex")]
    
    # the excel table contains two final tail rows which do not contain
    # tabular data, so we remove these
    variants = variants[1:(nrow(variants) - 2), ]
    
    # standardise the chrom, position and allele column names
    variants$chrom = gsub(" ", "", variants$Chromosome)
    variants$start_pos = gsub(" ", "", variants$Position..hg19.)
    variants$end_pos = variants$start_pos
    variants$ref_allele = gsub(" ", "", variants$Ref)
    variants$alt_allele = gsub(" ", "", variants$Allele)
    
    # sort out the "complex" allele events
    variants$alt_allele[61] = "S"
    variants$alt_allele[78] = "1D, -G"
    variants$alt_allele[116] = "1I, +C"
    variants$alt_allele[133] = "R"
    variants$alt_allele[185] = "1D, -G"
    
    # fix the alleles and positions for insertions and deletions
    deletions = grep("D, *-", variants$alt_allele)
    insertions = grep("I, *\\+", variants$alt_allele)
    
    # find the reference sequence at the site. Deletions use this as the
    # alternate allele, whereas the insertions use this as the reference allele
    alt_dels = apply(variants[deletions, ], 1, get_sequence_in_region)
    ref_ins = apply(variants[insertions, ], 1, get_sequence_in_region)
    
    # find the sequence at the site + the distance of the deletion
    variants$ref_allele[deletions] = paste(alt_dels, sapply(strsplit(variants$alt_allele[deletions], "D, *-"), "[", 2), sep="")
    variants$alt_allele[deletions] = alt_dels
    variants$ref_allele[insertions] = ref_ins
    variants$alt_allele[insertions] = paste(ref_ins, sapply(strsplit(variants$alt_allele[insertions], "I, *\\+"), "[", 2), sep="")
    
    # get the end coordinate, including those for the insertions and deletions
    variants$end_pos[deletions] = as.numeric(variants$end_pos[deletions]) + nchar(variants$ref_allele[deletions]) - 1
    
    variants = fix_het_alleles(variants)
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Person
    variants$study_code = "oroak_nature_2012"
    variants$publication_doi = "10.1038/nature10989"
    variants$study_phenotype = "autism"
    
    variants = merge(variants, gender, by="person_id", all.x=TRUE)
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}

#' get de novo data from the 2012 Iossifov et al autism exome study in Neuron
#'
#' Supplementary table 1 (where the non-coding SNVs have been excluded) and
#' supplementary table 2 from:
#' Iossifov et al. (2012) Neuron 74:285-299
#' doi: 10.1016/j.neuron.2012.04.009
#'
#' @return data frame of de novos, with standardised genome coordinates and VEP
#      consequences for each variant
iossifov_neuron_de_novos <- function() {
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0896627312003406/1-s2.0-S0896627312003406-mmc2.xlsx/272195/FULL/S0896627312003406/26c5ba3b72a2410ef43fec52a40f35e6/mmc2.xlsx"
    snvs = gdata::read.xls(url, sheet="SNV.v4.1-normlized", stringsAsFactors=FALSE)
    
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0896627312003406/1-s2.0-S0896627312003406-mmc4.xlsx/272195/FULL/S0896627312003406/6caa42b35609c2ed5910b5381ddd5335/mmc4.xlsx"
    indels = gdata::read.xls(url, sheet="ID.v4.1-normlized", stringsAsFactors=FALSE)
    
    # trim out the low quality de novos (as defined by a flag in the table)
    snvs = snvs[snvs$SNVFilter == 1, ]
    indels = indels[indels$IndelFilter == 1, ]
    
    # merge the SNV and indel de novo calls
    snvs = subset(snvs, select = c("quadId", "location", "variant", "effectGenes", "effectType", "inChild"))
    indels = subset(indels, select = c("quadId", "location", "variant", "effectGenes", "effectType", "inChild"))
    variants = rbind(snvs, indels)
    
    # get the coordinates and VEP consequence
    variants = fix_coordinates_with_allele(variants, "location", "variant")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    # exclude de novos not located within the coding sequence of a gene
    variants = variants[!(variants$consequence %in% NONCODING_CONSEQUENCES), ]
    
    variants$person_id = variants$quadId
    variants$study_code = "iossifov_neuron_2012"
    variants$publication_doi = "10.1016/j.neuron.2012.04.009"
    variants$study_phenotype = "autism"
    
    variants$sex = substr(variants$inChild, 4, 4)
    variants$sex[variants$sex == "M"] = "male"
    variants$sex[variants$sex == "F"] = "female"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}

#' get de novo data from the 2014 Iossifov et al. autism exome study in Nature
#'
#' De novo mutation data sourced from Supplementary table 2:
#' Iossifov et al. (2014) Nature 498:216-221
#' doi: 10.1038/nature13908
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
iossifov_nature_de_novos <- function() {
    tmpdir = tempdir()
    path = tempfile(tmpdir=tmpdir)
    
    # obtain the dataframe of de novo variants
    download.file("http://www.nature.com/nature/journal/v515/n7526/extref/nature13908-s2.zip", path)
    unzip(path, files="nature13908-s2/Supplementary Table 2.xlsx", exdir=tmpdir)
    variants = gdata::read.xls(file.path(tmpdir, "nature13908-s2", "Supplementary Table 2.xlsx"), stringsAsFactors=FALSE)
    unlink(path)
    
    variants = fix_coordinates(variants, "location", "vcfVariant")
    
    # exclude the de novos from the unaffected sibs
    variants = variants[!grepl("^s", variants$inChild), ]
    
    # get the sex of the probands
    variants$sex = substr(variants$inChild, 2, 2)
    variants$sex[variants$sex == "M"] = "male"
    variants$sex[variants$sex == "F"] = "female"
    
    # NOTE: the variant with the most severe consequence might not necessarily
    # NOTE: be within the listed gene. I haven't accounted for this yet.
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    variants = variants[!(variants$consequence %in% NONCODING_CONSEQUENCES), ]
    
    variants$person_id = variants$familyId
    variants$study_code = "iossifov_nature_2014"
    variants$publication_doi = "10.1038/nature13908"
    variants$study_phenotype = "autism"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}

#' get de novo data from the 2014 De Rubeis et al. autism exome study in Nature
#'
#' De novo mutation data sourced from Supplementary table 3:
#' De Rubeis et al. (2013) Nature 515:209-215
#' doi: 10.1038/nature13772
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
de_rubeis_de_novos <- function() {
    url = "http://www.nature.com/nature/journal/v515/n7526/extref/nature13772-s4.xlsx"
    
    variants = gdata::read.xls(url, sheet="De Novo", stringsAsFactors=FALSE)
    
    # exclude the final row, which contains a footnote
    variants = variants[1:(nrow(variants) - 1), ]
    
    # rename columns to match the other de novo datasets, and strip whitespace
    variants$start_pos = gsub("[ \t]", "", variants$Pos)
    variants$chrom = gsub("[ \t]", "", variants$Chr)
    variants$person_id = gsub("[ \t]", "", variants$Child_ID)
    variants$ref_allele = gsub("[ \t]", "", variants$Ref)
    variants$alt_allele = gsub("[ \t]", "", variants$Alt)
    
    # get the end position
    variants$end_pos = as.character(as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1)
    
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    # figure out the sex for each variant
    variants$sex = variants$Child_Sex
    variants$sex[variants$sex == "2"] = "female"
    variants$sex[variants$sex == "1"] = "male"
    
    variants$study_code = "derubeis_nature_2014"
    variants$publication_doi = "10.1038/nature13772"
    variants$study_phenotype = "autism"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
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
#' Supplementary table 1 (where the non-coding SNVs have been excluded) and
#' supplementary table 2 from:
#' Iossifov et al. (2012) Neuron 74:285-299
#' doi: 10.1016/j.neuron.2012.04.009
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
autism_de_novos <- function() {
    
    sanders = sanders_de_novos()
    oroak = oroak_de_novos()
    iossifov_neuron = iossifov_neuron_de_novos()
    iossifov_nature = iossifov_nature_de_novos()
    derubeis_nature = de_rubeis_de_novos()
    
    # exclude de novos identified in previous studies
    key_1 = paste(iossifov_neuron$person_id, iossifov_neuron$start_pos)
    key_2 = paste(iossifov_nature$person_id, iossifov_nature$start_pos)
    iossifov_nature = iossifov_nature[!(key_2 %in% key_1), ]
    
    autism_de_novos = rbind(sanders, oroak, iossifov_neuron, iossifov_nature,
        derubeis_nature)
    
    # remove de novos that have been dupicated between studies. These are easily
    # spotted as individuals who have IDs that are nearly identical between
    # different studies eg 14323.p1 vs 14323
    temp = autism_de_novos
    temp$person_id = sapply(strsplit(temp$person_id, "\\."), "[", 1)
    to_remove =  temp[duplicated(temp[, c("person_id", "chrom", "start_pos")]), ]
    autism_de_novos = autism_de_novos[!row.names(autism_de_novos) %in% row.names(to_remove), ]
    
    return(autism_de_novos)
}

#' get de novo data from Fromer et al. schizophrenia exome study
#'
#' De novo mutation data sourced from Supplementary table 1:
#' Fromer et al. (2014) Nature 506:179-184
#' doi: 10.1038/nature12929
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
fromer_de_novos <- function() {
    
    url = "http://www.nature.com/nature/journal/v506/n7487/extref/nature12929-s2.xlsx"
    variants = gdata::read.xls(url, stringsAsFactors=FALSE)
    
    variants$temp = strsplit(variants$Locus, ":")
    variants$chrom = sapply(variants$temp, "[", 1)
    variants$chrom = gsub("chr", "", variants$chrom)
    variants$start_pos = sapply(variants$temp, "[", 2)
    variants$end_pos = variants$start_pos
    
    # fix indel ranges
    indels = grepl("\\.\\.", variants$start_pos)
    variants$start_pos[indels] = sapply(strsplit(variants$start_pos[indels], "\\.\\."), "[", 1)
    variants$end_pos[indels] = sapply(strsplit(variants$end_pos[indels], "\\.\\."), "[", 2)
    
    variants$ref_allele = variants$Reference.allele
    variants$alt_allele = variants$Alternate.allele
    
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Proband.ID
    variants$study_code = "fromer_nature_2014"
    variants$publication_doi = "10.1038/nature12929"
    variants$study_phenotype = "schizophrenia"
    
    # get the sex for each variant
    variants$sex = variants$Proband.gender
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}

#' get de novo data from Zaidi et al. congenital heart disease exome study
#'
#' De novo mutation data sourced from Supplementary table 4:
#' Zaidi et al. (2013) Nature 498:220-223
#' doi: 10.1038/nature12141
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
zaidi_de_novos <- function() {
    
    url = "http://www.nature.com/nature/journal/v498/n7453/extref/nature12141-s1.pdf"
    
    # obtain the supplementary material
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    table_s4 = test$content[314:912]
    table_s4 = table_s4[table_s4 != ""]
    
    # clean up table S4
    table_s4 = gsub("^[ \t\f]+", "", table_s4) # trim leading whitespace
    table_s4 = gsub(" bp beyond exon ", "_bp_beyond_exon_", table_s4)
    table_s4 = gsub(" bp up of exon ", "_bp_up_of_exon_", table_s4)
    
    # drop the section breaks from the table, as well as the lines that are part
    # of the line breaks
    table_s4 = table_s4[!grepl("[Mm]utations", table_s4)]
    table_s4 = table_s4[!grepl("N A T U R E", table_s4)]
    table_s4 = table_s4[!grepl("SUPPLEMENTARY INFORMATION", table_s4)]
    
    # cull it down to the first few entries in each line
    split_strings = strsplit(table_s4, "[ \t]+")
    variants = data.frame(t(sapply(split_strings, "[", 1:11)))
    names(variants) = c("person_id", "category", "hgnc", "type", "aa_change", "dbSNP", "transcript", "protein", "chrom", "start_pos", "alleles")
    
    # drop out the controls
    variants = variants[variants$category != "Control", ]
    
    variants = fix_zaiidi_coordinates(variants, "alleles")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "zaidi_nature_2013"
    variants$publication_doi = "10.1038/nature12141"
    variants$study_phenotype = "congenital_heart_disease"
    variants$sex = NA
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}

rauch = rauch_de_novos()
deligt = deligt_de_novos()
deligt_gender = deligt$gender
deligt = deligt$variants
gilissen = gilissen_de_novos(deligt, deligt_gender)
epi4k = epi4k_de_novos()
autism = autism_de_novos()
fromer = fromer_de_novos()
zaidi = zaidi_de_novos()

published_de_novos = rbind(rauch, deligt, gilissen, epi4k, autism, fromer, zaidi)
published_de_novos$type = "snv"
published_de_novos$type[nchar(published_de_novos$ref_allele) != 1 | nchar(published_de_novos$alt_allele) != 1] = "indel"

# make sure all the columns are character
published_de_novos[] = lapply(published_de_novos, as.character)

save(published_de_novos, file="data/published_de_novos.rda", compress="xz")
# check that the gene names are fine

# Check that all the frameshifts are indels ()

# check that almost all the non-frameshifts are SNVs (aside from ones that have stop_gain consequences)
