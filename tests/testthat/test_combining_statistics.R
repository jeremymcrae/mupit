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
        hgnc  lof_p  func_p
        GENE1 0.01   0.01
        GENE2 0.0001 0.001 ")
    
    clust = read.table(header=TRUE, text="
        gene_id mutation_category probability
        GENE1   missense          0.1
        GENE1   nonsense          0.1
        GENE2   missense          0.001
        GENE2   nonsense          0.001")
    
    output = read.table(header=TRUE, text="
    hgnc  lof_p func_p p_missense_clust p_nonsense_clust p_combined    combined_fdr
    GENE1 1e-02 0.01   0.100            0.100            7.9077553e-03 1.000000000
    GENE2 1e-04 0.001  0.001            0.001            1.7118110e-06 0.031668477")
    
    expect_equal(combine_enrichment_and_clustering(enriched, clust), output)
    
    # and check that the values are correct when we alter the num_tests variable
    output = read.table(header=TRUE, text="
    hgnc  lof_p func_p p_missense_clust p_nonsense_clust p_combined   combined_fdr
    GENE1 1e-02 0.01   0.100            0.100            7.9077553e-03 0.015815510557
    GENE2 1e-04 0.001  0.001            0.001            1.7118110e-06 6.84723826e-06")
    
    expect_equal(combine_enrichment_and_clustering(enriched, clust, num_tests=4), output)
})
