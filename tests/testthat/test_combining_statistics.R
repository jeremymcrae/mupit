# unit testing for the mupit functions

library(mupit)
library(testthat)
library(jsonlite)

context("combining analyses checks")

test_that("fishersMethod is correct for simple vector", {
    p_values = c(0.01, 0.001)
    expect_equal(fishersMethod(p_values), 0.000125129254)
})

test_that("fishersMethod is correct for single value", {
    p_values = c(0.01)
    expect_equal(fishersMethod(p_values), 0.01)
})

test_that("fishersMethod is correct for vector with NA", {
    p_values = c(0.01, NA)
    expect_equal(fishersMethod(p_values), 0.01)
})

test_that("fishersMethod is correct for vector with NA", {
    p_values = c(0.01, 0.001, NA)
    expect_equal(fishersMethod(p_values), 0.000125129254)
})

test_that("fishersMethod is correct without any p-values", {
    p_values = c(NA, NA)
    expect_equal(fishersMethod(p_values), NA)
})

test_that("fishersMethod is correct for vector with zero", {
    p_values = c(0, 0.01)
    expect_equal(fishersMethod(p_values), 0)
})

test_that("combine_enrichment_and_clustering is correct for small table", {
    enriched = read.table(header=TRUE, text="
        hgnc  p_lof  p_func
        GENE1 0.01   0.01
        GENE2 0.0001 0.001")
    
    clust = read.table(header=TRUE, text="
        gene_id mutation_category probability
        GENE1   missense          0.1
        GENE1   nonsense          0.1
        GENE2   missense          0.001
        GENE2   nonsense          0.001")
    
    output = read.table(header=TRUE, text="
    hgnc  p_lof p_func p_missense_clust p_nonsense_clust p_combined        combined_fdr p_min
    GENE1 1e-02 0.01   0.100            0.100            7.90775527898e-03 1.000000000  0.00790775527898214
    GENE2 1e-04 0.001  0.001            0.001            1.48155105579e-05 0.274086945  1.48155105579643e-05 ")
    
    expect_equal(combine_enrichment_and_clustering(enriched, clust), output)
    
    # and check that the values are correct when we alter the num_tests variable
    output = read.table(header=TRUE, text="
    hgnc  p_lof p_func p_missense_clust p_nonsense_clust p_combined   combined_fdr         p_min
    GENE1 1e-02 0.010  0.100            0.100            7.90775527898e-03 0.0158155105579 0.00790775527898214
    GENE2 1e-04 0.001  0.001            0.001            1.48155105579e-05 5.926204223e-05 1.48155105579643e-05")
    
    expect_equal(combine_enrichment_and_clustering(enriched, clust, num_tests=4), output)
})
