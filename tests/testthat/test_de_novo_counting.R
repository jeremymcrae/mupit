# unit testing for the mupit functions

library(mupit)
library(testthat)

context("De novo counting checks")

test_that("get_de_novo_counts output is correct for simplest table", {
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_1 ARID1B 6 157431695 missense_variant snv
        person_2 ARID1B 6 157502190 stop_gained snv")
    
    output = read.table(header=TRUE, stringsAsFactors=FALSE, text="
        hgnc   lof_snv missense_snv chrom min_pos   lof_indel missense_indel
        ARID1B 1       1            6     157431695 0         0",
        colClasses=c("character", "numeric", "numeric", "character", "numeric",
            "numeric", "numeric"))
    
    expect_equal(get_de_novo_counts(vars), output)
})

test_that("get_de_novo_counts output is correct for a multi-gene dataset", {
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_1 ARID1B 6 157150547 inframe_deletion indel
        person_2 ARID1B 6 157431695 missense_variant snv
        person_3 ARID1B 6 157454186 frameshift_variant indel
        person_4 ARID1B 6 157502190 stop_gained snv
        person_4 KMT2A 11 118367083 splice_donor_variant snv")
    
    output = read.table(header=TRUE, text="
        hgnc   lof_indel lof_snv missense_indel missense_snv chrom min_pos
        ARID1B 1       1            1         1              6     157150547
        KMT2A  0       1            0         0              11    118367083",
        colClasses=c("character", "numeric", "numeric", "numeric", "numeric",
            "character", "numeric"))
    
    expect_equal(get_de_novo_counts(vars), output)
})

test_that("get_de_novo_counts output is correct without functional variants", {
    # check synonymous_variant gained on its own raises an error
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 synonymous_variant snv")
    expect_error(get_de_novo_counts(vars), "replacement has 0 rows, data has 1")
    
    # check that non-functional, or variants not included in our model don't
    # contribute to the gene counts
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_1 ARID1B 6 157431695 missense_variant snv
        person_2 ARID1B 6 157502190 synonymous_variant snv
        person_2 ARID1B 6 157502190 splice_region_variant snv
        person_2 ARID1B 6 157502190 intron_variant snv
        ")
    
    output = read.table(header=TRUE, stringsAsFactors=FALSE, text="
        hgnc   missense_snv chrom min_pos   lof_indel missense_indel lof_snv
        ARID1B 1            6     157431695 0         0              0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    
    expect_equal(get_de_novo_counts(vars), output)
})

test_that("get_de_novo_counts output is correct for LoF consequence types", {
    # check stop_gained
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 stop_gained snv")
    
    output = read.table(header=TRUE, text="
        hgnc   lof_snv chrom min_pos   lof_indel missense_indel missense_snv
        ARID1B 1       6     157150547 0         0       0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check stop_gained")
    
    # check splice_acceptor_variant
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 splice_acceptor_variant snv")
    
    output = read.table(header=TRUE, text="
        hgnc   lof_snv chrom min_pos   lof_indel missense_indel missense_snv
        ARID1B 1       6     157150547 0         0       0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check splice_acceptor_variant")
    
    # check splice_donor_variant
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 splice_donor_variant snv")
    
    output = read.table(header=TRUE, text="
        hgnc   lof_snv chrom min_pos   lof_indel missense_indel missense_snv
        ARID1B 1       6     157150547 0         0       0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check splice_donor_variant")
    
    # check frameshift_variant
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 frameshift_variant indel")
    
    output = read.table(header=TRUE, text="
        hgnc   lof_indel chrom min_pos   missense_indel lof_snv missense_snv
        ARID1B 1         6     157150547 0              0       0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check frameshift_variant")
})

test_that("get_de_novo_counts output is correct for functional consequence types", {
    # check missense_variant
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 missense_variant snv")
    
    output = read.table(header=TRUE, text="
        hgnc   missense_snv chrom min_pos   lof_indel missense_indel lof_snv
        ARID1B 1            6     157150547 0         0              0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check missense_variant")
    
    # check initiator_codon_variant
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 initiator_codon_variant snv")
    
    output = read.table(header=TRUE, text="
        hgnc   missense_snv chrom min_pos   lof_indel missense_indel lof_snv
        ARID1B 1            6     157150547 0         0              0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check initiator_codon_variant")
    
    # check stop_lost
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 stop_lost snv")
    
    output = read.table(header=TRUE, text="
        hgnc   missense_snv chrom min_pos   lof_indel missense_indel lof_snv
        ARID1B 1            6     157150547 0         0              0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check stop_lost")
    
    # check coding_sequence_variant
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 coding_sequence_variant snv")
    
    output = read.table(header=TRUE, text="
        hgnc   missense_snv chrom min_pos   lof_indel missense_indel lof_snv
        ARID1B 1            6     157150547 0         0              0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check coding_sequence_variant")
    
    # check inframe_deletion
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 inframe_deletion indel")
    
    output = read.table(header=TRUE, text="
        hgnc   missense_indel chrom min_pos   lof_indel lof_snv missense_snv
        ARID1B 1              6     157150547 0         0       0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check inframe_deletion")
    
    # check inframe_insertion
    vars = read.table(header=TRUE, colClasses="character", text="
        person_id hgnc chrom start_pos consequence type
        person_3 ARID1B 6 157150547 inframe_insertion indel")
    
    output = read.table(header=TRUE, text="
        hgnc   missense_indel chrom min_pos   lof_indel lof_snv missense_snv
        ARID1B 1              6     157150547 0         0       0",
        colClasses=c("character", "numeric", "character", "numeric", "numeric",
            "numeric", "numeric"))
    expect_equal(get_de_novo_counts(vars), output, label="check inframe_insertion")
})
