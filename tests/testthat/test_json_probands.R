# unit testing for the mupit functions

library(mupit)
library(testthat)
library(jsonlite)

context("JSON output checks")

test_that("write_probands_by_gene output is correct", {
    
    de_novos = read.table(header=TRUE, text="
        person_id  hgnc
        person_1   GENE1
        person_1   GENE2
        person_2   GENE1",
        colClasses=c("character", "character"))
    
    # write the file
    path = tempfile()
    write_probands_by_gene(de_novos, path)
    
    # check that the function output matches our expectations
    json = jsonlite::fromJSON(path)
    expect_true(file.exists(path))
    expect_equal(json, list("GENE1"=c("person_1", "person_2"), "GENE2"=c("person_1")))
    
    # delete the temporary output file
    unlink(path)
})
