# unit testing for the mupit functions

library(mupit)
library(testthat)

context("Loading DDD de novos checks")

test_that("standardise_ddd_de_novos output is correct", {
    
    # define a text file for loading, and write it to a file
    text = c("person_stable_id\tsex\tchrom\tpos\tref\talt\tsymbol\tconsequence\tvar_type",
        "person_1\tF\t1\t1000\tA\tG\tGENE1\tmissense_variant\tDENOVO-SNP",
        "person_2\tM\t1\t1010\tC\tGGG\tGENE2\tframeshift_variant\tDENOVO-INDEL",
        "person_3\tF\tX\t2000\tCCC\tG\tGENE3\tframeshift_variant\tDENOVO-INDEL")
    path = tempfile()
    writeLines(text, path, sep="\n")
    
    # define the expected output
    output = read.table(header=TRUE, text="
        person_id sex chrom start_pos end_pos ref_allele alt_allele hgnc  consequence       study_code  publication_doi study_phenotype type
        person_1  F   1     1000      1000    A          G          GENE1 missense_variant  ddd_unpublished  NA  developmental_disorders snv
        person_2  M   1     1010      1010    C          GGG        GENE2 frameshift_variant ddd_unpublished  NA  developmental_disorders indel
        person_3  F   X     2000      2002    CCC        G          GENE3 frameshift_variant ddd_unpublished  NA  developmental_disorders indel",
        colClasses=c("character", "character", "character", "numeric",
            "numeric", "character", "character", "character", "character",
            "character", "logical", "character", "character"))
    
    expect_equal(standardise_ddd_de_novos(path), output)
    unlink(path)
})

test_that("get_ddd_de_novos output is correct", {
    
    # define a text file for loading, and write it to a file
    text = c("person_stable_id\tsex\tchrom\tpos\tref\talt\tsymbol\tconsequence\tvar_type",
        "person_1\tF\t1\t1000\tA\tG\tGENE1\tmissense_variant\tDENOVO-SNP",
        "person_2\tM\t1\t1010\tC\tGGG\tGENE2\tframeshift_variant\tDENOVO-INDEL",
        "person_3\tF\tX\t2000\tCCC\tG\tGENE3\tframeshift_variant\tDENOVO-INDEL")
    path = tempfile()
    writeLines(text, path, sep="\n")
    
    # define the expected output
    output = read.table(header=TRUE, text="
        person_id sex chrom start_pos end_pos ref_allele alt_allele hgnc  consequence       study_code  publication_doi study_phenotype type
        person_1  F   1     1000      1000    A          G          GENE1 missense_variant  ddd_unpublished  NA  developmental_disorders snv
        person_2  M   1     1010      1010    C          GGG        GENE2 frameshift_variant ddd_unpublished  NA  developmental_disorders indel
        person_3  F   X     2000      2002    CCC        G          GENE3 frameshift_variant ddd_unpublished  NA  developmental_disorders indel",
        colClasses=c("character", "character", "character", "numeric",
            "numeric", "character", "character", "character", "character",
            "character", "logical", "character", "character"))
    
    expect_equal(get_ddd_de_novos(path), output)
    unlink(path)
})

test_that("get_ddd_de_novos output is correct when removing the diagnosed set", {
    
    # define a text file for loading, and write it to a file
    text = c("person_stable_id\tsex\tchrom\tpos\tref\talt\tsymbol\tconsequence\tvar_type",
        "person_1\tF\t1\t1000\tA\tG\tGENE1\tmissense_variant\tDENOVO-SNP",
        "person_2\tM\t1\t1010\tC\tGGG\tGENE2\tframeshift_variant\tDENOVO-INDEL",
        "person_3\tF\tX\t2000\tCCC\tG\tGENE3\tframeshift_variant\tDENOVO-INDEL")
    path = tempfile()
    writeLines(text, path, sep="\n")
    
    exclude = data.frame("person_id"=c("person_1", "person_2"), "sex"=c("F", "M"))
    
    # define the expected output
    output = read.table(header=TRUE, text="
        person_id sex chrom start_pos end_pos ref_allele alt_allele hgnc  consequence       study_code  publication_doi study_phenotype type
        person_3  F   X     2000      2002    CCC        G          GENE3 frameshift_variant ddd_unpublished  NA  developmental_disorders indel",
        colClasses=c("character", "character", "character", "numeric",
            "numeric", "character", "character", "character", "character",
            "character", "logical", "character", "character"))
    
    # set the row names for the output, taking into account the saples that have
    # been removed.
    attr(output, "row.names") = as.integer(3)
    
    expect_equal(get_ddd_de_novos(path, diagnosed=exclude), output)
    unlink(path)
})

test_that("get_ddd_de_novos output is correct when selecting a subset", {
    
    # define a text file for loading, and write it to a file
    text = c("person_stable_id\tsex\tchrom\tpos\tref\talt\tsymbol\tconsequence\tvar_type",
        "person_1\tF\t1\t1000\tA\tG\tGENE1\tmissense_variant\tDENOVO-SNP",
        "person_2\tM\t1\t1010\tC\tGGG\tGENE2\tframeshift_variant\tDENOVO-INDEL",
        "person_3\tF\tX\t2000\tCCC\tG\tGENE3\tframeshift_variant\tDENOVO-INDEL")
    path = tempfile()
    writeLines(text, path, sep="\n")
    
    include = c("person_1")
    
    # define the expected output
    output = read.table(header=TRUE, text="
        person_id sex chrom start_pos end_pos ref_allele alt_allele hgnc  consequence       study_code  publication_doi study_phenotype type
        person_1  F   1     1000      1000    A          G          GENE1 missense_variant  ddd_unpublished  NA  developmental_disorders snv",
        colClasses=c("character", "character", "character", "numeric",
            "numeric", "character", "character", "character", "character",
            "character", "logical", "character", "character"))
    
    expect_equal(get_ddd_de_novos(path, subset=include), output)
    unlink(path)
})

test_that("load_ddg2p output is correct", {
    
    text = c("chr\tstart\tstop\tgene\ttype\tmode\tmech",
        "1\t100\t150\tGENE1\tProbable DD gene\tBiallelic\tLoss of function\tSyndromeZ",
        "1\t200\t250\tGENE2\tConfirmed DD Gene\tBiallelic\tLoss of function\tSyndromeZ",
        "1\t300\t350\tGENE3\tPossible DD Gene\tMonoallelic\tLoss of function\tSyndromeZ",
        "1\t400\t450\tGENE4\tBoth DD and IF\tMonoallelic\tLoss of function\tSyndromeZ",
        "1\t500\t550\tGENE5\tBoth DD and IF\tBialleic\tLoss of function\tSyndromeZ",
        "X\t600\t650\tGENE6\tBoth DD and IF\tX-linked dominant\tLoss of function\tSyndromeZ",
        "X\t700\t750\tGENE7\tBoth DD and IF\tHemizygous\tLoss of function\tSyndromeZ")
    
    path = tempfile()
    writeLines(text, path, sep="\n")
    
    # define the expected output
    output = read.table(sep="\t", strip.white=TRUE, header=TRUE, text="
        chr\tstart\tstop\tgene \ttype             \tmode             \tmech            \tdescription\tdominant\themizygous
        1  \t100  \t150 \tGENE1\tProbable DD gene \tBiallelic        \tLoss of function\tSyndromeZ  \tFALSE   \tFALSE
        1  \t200  \t250 \tGENE2\tConfirmed DD Gene\tBiallelic        \tLoss of function\tSyndromeZ  \tFALSE   \tFALSE
        1  \t400  \t450 \tGENE4\tBoth DD and IF   \tMonoallelic      \tLoss of function\tSyndromeZ  \tTRUE    \tFALSE
        1  \t500  \t550 \tGENE5\tBoth DD and IF   \tBialleic         \tLoss of function\tSyndromeZ  \tFALSE   \tFALSE
        X  \t600  \t650 \tGENE6\tBoth DD and IF   \tX-linked dominant\tLoss of function\tSyndromeZ  \tTRUE    \tFALSE
        X  \t700  \t750 \tGENE7\tBoth DD and IF   \tHemizygous       \tLoss of function\tSyndromeZ  \tFALSE   \tTRUE",
        colClasses=c("character", "numeric", "numeric", "character",
            "character", "character", "character", "character", "logical", "logical"))
    
    expect_equal(load_ddg2p(path), output)
    unlink(path)
})
