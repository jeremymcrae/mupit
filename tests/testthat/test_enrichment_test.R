# unit testing for the mupit functions

library(mupit)
library(testthat)

context("Enrichment testing checks")

test_that("test_enrichment output is correct", {
    counts = read.table(header=TRUE, text="
        hgnc   lof_indel lof_snv missense_indel missense_snv chrom min_pos
        ARID1B 1       1            1         1              6     157150547
        KMT2A  0       0            0         1              11    118367083")
    
    rates = read.table(header=TRUE, text="
        hgnc   chrom snv.lof.rate indel.lof.rate snv.missense.rate indel.missense.rate
        ARID1B 6     0.01        0.005          0.05             0.005
        KMT2A  11    0.01        0.005          0.05             0.005")
    
    # define the expected output, including the exact P-values expected from the
    # submitted data frame inputs
    output = read.table(header=TRUE, text="
        hgnc   chrom lof_indel lof_snv missense_indel missense_snv min_pos   p_lof        p_func
        ARID1B 6     1         1       1              1            157150547 0.000111381302 9.459951611e-07
        KMT2A  11    0         0       0              1            118367083 1.000000000000 6.760618009e-02")
    
    expect_equal(gene_enrichment(rates, counts), output)
    
    # and check that it can also pick up the other gene in the rate table
    counts = read.table(header=TRUE, text="
        hgnc   lof_indel lof_snv missense_indel missense_snv chrom min_pos
        KMT2A  0       0            0         1              11    118367083")
    output = read.table(header=TRUE, text="
        hgnc   chrom lof_indel lof_snv missense_indel missense_snv min_pos   p_lof p_func
        ARID1B 6     NA        NA      NA             NA           NA        NA    NA
        KMT2A  11    0         0       0              1            118367083 1     6.760618009e-02")
    expect_equal(gene_enrichment(counts, rates), output)
})

test_that("analyse_gene_enrichment output is correct", {
    trios = list(male=1000, female=1000)
    vars = read.table(header=TRUE, text="
        person_id hgnc chrom start_pos consequence type
        person_1 ARID1B 6 157150547 inframe_deletion indel
        person_2 ARID1B 6 157431695 missense_variant snv
        person_3 ARID1B 6 157454186 frameshift_variant indel
        person_4 ARID1B 6 157502190 stop_gained snv
        person_4 KMT2A 11 118367083 missense_variant snv")
    
    rates = read.table(header=TRUE, text="
        hgnc    chrom  syn  mis  non   splice_site  frameshift
        ARID1B  6      -5   -6   -6.5  -5           -7
        KMT2A   11     -5   -6   -6.5  -5           -7")
    
    # define the expected output, including the exact P-values expected from the
    # submitted data frame inputs
    output = read.table(header=TRUE, text="
        hgnc   chrom lof_indel lof_snv missense_indel missense_snv p_lof         p_func
        ARID1B 6     1         1       1              1            0.0008401753  1.7364753290e-07
        KMT2A  11    0         0       0              1            1.0           4.4572185724e-02",
        colClasses=c("character", "numeric", "numeric", "numeric", "numeric",
            "numeric", "numeric", "numeric"))
    
    expect_equal(analyse_gene_enrichment(vars, trios, rates=rates), output)
    
    # and check that we can generate a Manhattan-style plot. This just checks
    # that we don't get any errors, it doesn't check if the plot is correct
    pdf_path = tempfile(fileext=".pdf")
    out = analyse_gene_enrichment(vars, trios, plot_path=pdf_path, rates=rates)
    expect_true(file.exists(pdf_path))
    unlink(pdf_path)
})
