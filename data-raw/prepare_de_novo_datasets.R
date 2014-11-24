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
    table_s2_text = test$content[769:847]
    table_s3_text = test$content[879:890]
    
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
    
    variants$study_code = "rauch_lancet_2012"
    variants$publication_doi = "10.1016/S0140-6736(12)61480-9"
    variants$study_phenotype = "intellectual_disability"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
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
    # path = tempfile()
    path = "corrupted_pdf.pdf"
    system(paste("wget --referer=", url, " --cookies=on --load-cookies=", cookie, " --keep-session-cookies --save-cookies=", cookie, " ", url, " -O ", path, sep=""))
    # download.file(url, path, extra=c("--cookies=on", "--load-cookies=cookie.txt", "--keep-session-cookies", "--save-cookies=cookie.txt"))
    
    # repair the pdf with ghostscript
    # repaired = tempfile()
    repaired = "repaired.pdf"
    system(paste("gs -o ", repaired, " -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress ", path, sep=""))
    
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri=repaired), language = "en", id = "id1")
    
    # delet the pdf and temporary files
    unlink(temp)
    unlink(cookie)
    unlink(path)
    unlink(repaired)
    
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
    # variants$consequence = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "deligt_nejm_2012"
    variants$publication_doi = "10.1056/NEJMoa1206524"
    variants$study_phenotype = "intellectual_disability"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence", 
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
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
gilissen_de_novos <- function() {
    
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
    
    # split the tabel, and drop the erroneous lines from the bad line
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
    # variants$consequence = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "gilissen_nature_2014"
    variants$publication_doi = "10.1038/nature13394"
    variants$study_phenotype = "intellectual_disability"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence", 
        "study_code", "publication_doi", "study_phenotype"))
    
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
    
    variants = gdata::read.xls("http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0002929714003838/1-s2.0-S0002929714003838-mmc2.xlsx/276895/FULL/S0002929714003838/bf21945d72e3297fc44969dc0296f4f1/mmc2.xlsx", stringsAsFactors=FALSE)
    
    # the excel table contains three final tail rows which do not contain
    # tabular data, so we remove these
    variants = variants[1:(nrow(variants) - 3), ]
    
    variants = fix_coordinates_with_allele(variants, 
        "hg19.coordinates..chr.position.", "Ref.Alt.alleles")
    
    # get the HGNC symbol
    # NOTE: the variant with the most severe consequence might not necessarily  
    # NOTE: be within the listed gene. I haven't accounted for this yet.
    # variants$consequence = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    
    # get the hgnc symbol, and clean any anomalies
    # variants$hgnc = variants$Gene
    # variants$hgnc = gsub(" \\(MLL4\\)", "", variants$hgnc)
    # variants$hgnc = gsub(" \\^", "", variants$hgnc)
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Child.ID
    variants$study_code = "epi4k_ajhg_2014"
    variants$publication_doi = "10.1016/j.ajhg.2014.08.013"
    variants$study_phenotype = "epilepsy"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
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
    variants$hgnc = NA
    for (row_num in 1:nrow(variants)) {
        # variants$consequence[row_num] = get_vep_consequence(variants[row_num, ], verbose=TRUE)
        vep = get_vep_consequence(variants[row_num, ], verbose=TRUE)
        variants$consequence[row_num] = vep$consequence
        variants$hgnc[row_num] = vep$gene
    }
    # variants$hgnc = variants$Gene
    variants$person_id = variants$Child_ID
    variants$study_code = "sanders_nature_2012"
    variants$publication_doi = "10.1038/nature10945"
    variants$study_phenotype = "autism"
    
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
    snvs = subset(snvs, select = c("quadId", "location", "variant", "effectGenes", "effectType"))
    indels = subset(indels, select = c("quadId", "location", "variant", "effectGenes", "effectType"))
    variants = rbind(snvs, indels)
    
    # get the coordinates and VEP consequence
    variants = fix_coordinates_with_allele(variants, "location", "variant")
    # variants$consequence = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    # get hgnc symbols for the few genes that have missing values
    # no_gene = variants$effectGenes == ""
    # variants$effectGenes[no_gene] = apply(variants[no_gene, ], 1, get_gene_id_for_variant)
    
    # extract the hgnc symbol from a gene:consequence string
    # variants$hgnc = sapply(strsplit(variants$effectGenes, ":"), "[[", 1)
    
    # exclude de novos not located within the coding sequence of a gene
    variants = variants[!(variants$consequence %in% NONCODING_CONSEQUENCES), ]
    
    variants$person_id = variants$quadId
    variants$study_code = "iossifov_neuron_2012"
    variants$publication_doi = "10.1016/j.neuron.2012.04.009"
    variants$study_phenotype = "autism"
    
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
    
    # NOTE: the variant with the most severe consequence might not necessarily  
    # NOTE: be within the listed gene. I haven't accounted for this yet.
    # variants$consequence = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    variants = variants[!(variants$consequence %in% NONCODING_CONSEQUENCES), ]
    
    # get the HGNC symbol
    variants$person_id = variants$familyId
    variants$study_code = "iossifov_nature_2014"
    variants$publication_doi = "10.1038/nature13908"
    variants$study_phenotype = "autism"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
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
    
    # exclude the final row, which is contains a footnote
    variants = variants[1:(nrow(variants) - 1), ]
    
    # rename columns to match the other de novo datasets, and strip whitespace
    variants$start_pos = gsub("[ \t]", "", variants$Pos)
    variants$chrom = gsub("[ \t]", "", variants$Chr)
    variants$person_id = gsub("[ \t]", "", variants$Child_ID)
    variants$ref_allele = gsub("[ \t]", "", variants$Ref)
    variants$alt_allele = gsub("[ \t]", "", variants$Alt)
    
    # get the end position
    variants$end_pos = as.character(as.numeric(variants$start_pos) + nchar(variants$ref_allele))
    # variants$consequence = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "derubeis_nature_2014"
    variants$publication_doi = "10.1038/nature13772"
    variants$study_phenotype = "autism"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
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
    iossifov_nature = iossifov_nature[key_2 %in% key_1, ]
    
    # TODO: exclude derubeis_nature variants from previous studies
    
    autism_de_novos = rbind(sanders, oroak, iossifov_neuron, iossifov_nature,
        derubeis_nature)
    
    return(autism_de_novos)
    
    # autism_de_novos = read.delim(file.path("data-raw", "de_novo_datasets", "autism_v3_PJ.txt"), header=TRUE, colClasses = "character")
    
    # # select only de novos in probands
    # autism_de_novos = autism_de_novos[which(autism_de_novos$pheno == "Pro"), ]
    
    # # standardise the SNV or INDEL flag
    # TYPE.index = which(nchar(autism_de_novos$ref.1) != nchar(autism_de_novos$var))
    # autism_de_novos$TYPE = "INDEL"
    # autism_de_novos$TYPE[TYPE.index] = "SNV"
    
    # # standardise the columns, and column names
    # autism_de_novos = subset(autism_de_novos, select = c("INFO.HGNC", "INFO.CQ", "pos", "CHROM", "TYPE"))
    # names(autism_de_novos) = c("hgnc", "consequence", "position", "chrom", "type")
    
    # # TODO: what I actually want is: 
    # # person_id, chrom, start_pos, end_pos, ref_allele, alt_allele, hgnc_symbol, consequence, variant_type, study_code, publication_doi, and study_phenotype
    # autism_de_novos$study = "autism"
    
    # save(autism_de_novos, file="data/autism_de_novos.rda")
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
    
    # variants$hgnc = variants$Genes
    variants$person_id = variants$Proband.ID
    variants$study_code = "fromer_nature_2014"
    variants$publication_doi = "10.1038/nature12929"
    variants$study_phenotype = "schizophrenia"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
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
    
    variants$study_code = "zaiidi_nature_2013"
    variants$publication_doi = "10.1038/nature12141"
    variants$study_phenotype = "congenital_heart_disease"
    
    variants = subset(variants, select=c("person_id", "chrom", "start_pos", 
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence", 
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}

rauch = rauch_de_novos()
deligt = deligt_de_novos()
gilissen = gilissen_de_novos()
epi4k = epi4k_de_novos()
autism = autism_de_novos()
fromer = fromer_de_novos()
zaiidi = zaidi_de_novos()

published_de_novos = rbind(rauch, deligt, gilissen, epi4k, autism, fromer, zaiidi)
published_de_novos$type = "snv"
published_de_novos$type[nchar(published_de_novos$ref_allele) != 1 | nchar(published_de_novos$alt_allele) != 1] = "indel"

save(published_de_novos, file="../data/published_de_novos.rda", compress="xz")
# check that the gene names are fine

# Check that all the frameshifts are indels ()

# check that almost all the non-frameshifts are SNVs (aside from ones that have stop_gain consequences)

