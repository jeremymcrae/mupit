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
        ARID1B 6     1         1       1              1            157150547 0.0001108251 9.327823173e-07
        KMT2A  11    0         0       0              1            118367083 0.9851119396 6.526756739e-02")
    
    expect_equal(test_enrichment(rates, counts), output)
})
