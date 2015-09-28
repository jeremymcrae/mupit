# unit testing for the mupit functions

library(mupit)
library(testthat)

context("Mutation rate checks")

test_that("get_gene_based_mutation_rates output is correct", {
    trios = list(male=1000, female=1000)
    rates = read.table(header=TRUE, text="
        hgnc    chrom  syn  mis  non   splice_site  frameshift
        ARID1B  6      -5   -6   -6.5  -5           -7
        KMT2A   11     -5   -6   -6.5  -5           -7
        AMELX   X      -5   -6   -6.5  -5           -7")
    
    # define the expected output, including the exact rates expected from the
    # submitted data frame inputs
    output = read.table(header=TRUE, text="
        hgnc    chrom  snv.missense.rate  snv.lof.rate   indel.missense.rate  indel.lof.rate
        ARID1B  6      0.004000000000     0.04126491106  3.311546841e-05      0.0002980392157
        KMT2A   11     0.004000000000     0.04126491106  3.311546841e-05      0.0002980392157
        AMELX   X      0.002454545455     0.02532164997  2.032085561e-05      0.0001828877005")
    
    expect_equal(get_gene_based_mutation_rates(trios, rates), output)
})

test_that("get_gene_based_mutation_rates output is correct for the default gene rates", {
    trios = list(male=1000, female=1000)
    
    # define the expected output, including the exact rates expected from the
    # submitted data frame inputs
    output = read.table(header=TRUE, text="
        hgnc   chrom snv.missense.rate  snv.lof.rate    indel.missense.rate  indel.lof.rate
        A1BG   19    0.06229597556      0.003251732050  0.0003121266868      0.002809140181
        A1CF   10    0.06142709700      0.007483523281  0.0003794293330      0.003414863997
        A2M    12    0.15705797416      0.008708833575  0.0009281772135      0.008353594922
        A2ML1  12    0.15774659415      0.014581549019  0.0009156529573      0.008240876616
        A4GALT 22    0.06354186994      0.001624112908  0.0002227567651      0.002004810886
        A4GNT  3     0.03749384601      0.002070752842  0.0002146005854      0.001931405269",
        colClasses=c("factor", "character", "numeric", "numeric", "numeric",
            "numeric"))
    
    result = get_gene_based_mutation_rates(trios)
    # make sure that we get the correct number of rows
    expect_equal(nrow(result), 18271)
    
    # restrict the gene rates to the top few rows of the table, and only check
    # that those values come out exactly correct
    result = head(result)
    result$hgnc = factor(result$hgnc)
    expect_equal(result, output)
})

test_that("get_gene_based_mutation_rates output is correct for the default gene rates", {
    trios = list(male=1000, female=1000)
    
    # define the expected output, including the exact rates expected from the
    # submitted data frame inputs
    output = read.table(header=TRUE, text="
    hgnc     chrom snv.missense.rate snv.lof.rate indel.missense.rate indel.lof.rate
    TSPAN6   X     0.017925248       0.0013178331 9.6007090e-05       0.0008640638
    TNMD     X     0.023171662       0.0017035404 0.00012410672       0.0011169605
    DPM1     20    0.035069652       0.002578260  1.878320e-04        0.0016904880
    SCYL3    1     0.088228278       0.006486390  4.725480e-04        0.0042529320
    C1orf112 1     0.101409084       0.007455420  5.431440e-04        0.0048882960
    FGR      1     0.06293538        0.004626900  3.370800e-04        0.0030337200",
        colClasses=c("factor", "character", "numeric", "numeric", "numeric",
            "numeric"))
    
    result = get_length_based_rates(trios)
    # make sure that we get the correct number of rows
    expect_equal(nrow(result), 20823)
    
    # restrict the gene rates to the top few rows of the table, and only check
    # that those values come out exactly correct
    result = head(result)
    result$hgnc = factor(result$hgnc)
    expect_equal(result, output)
})

test_that("correct_for_x_chrom output is correct", {
    trios = list(male=1000, female=1000)
    rates = read.table(header=TRUE, text="
        hgnc    chrom  snv.missense.rate  snv.lof.rate  indel.missense.rate  indel.lof.rate
        ARID1B  6      0.001              0.001         0.001                0.001
        KMT2A   11     0.001              0.001         0.001                0.001
        AMELX   X      0.001              0.001         0.001                0.001")
    
    # define the expected output, including the exact rates expected from the
    # submitted data frame inputs
    output = read.table(header=TRUE, text="
        hgnc    chrom  snv.missense.rate  snv.lof.rate     indel.missense.rate  indel.lof.rate
        ARID1B  6      0.001              0.001            0.001                0.001
        KMT2A   11     0.001              0.001            0.001                0.001
        AMELX   X      0.0006136363636    0.0006136363636  0.0006136363636      0.0006136363636")
    
    expect_equal(correct_for_x_chrom(rates, trios$male, trios$female), output, label="balanced group of males and females")
    expect_equal(correct_for_x_chrom(rates, 0, trios$female), rates, label="chrX adjustment with no males")
})

test_that("adjust_indel_rates output is correct", {
    rates = read.table(header=TRUE, text="
        hgnc    chrom  snv.missense.rate  snv.lof.rate  indel.missense.rate  indel.lof.rate
        ARID1B  6      0.001              0.001         0.001                0.001
        KMT2A   11     0.001              0.001         0.001                0.001
        AMELX   X      0.001              0.001         0.001                0.001")
    
    # define the expected output, including the exact rates expected from the
    # submitted data frame inputs
    output = read.table(header=TRUE, text="
        hgnc    chrom  snv.missense.rate  snv.lof.rate     indel.missense.rate  indel.lof.rate
        ARID1B  6      0.001              0.001            0.0007450980392      0.0007450980392
        KMT2A   11     0.001              0.001            0.0007450980392      0.0007450980392
        AMELX   X      0.001              0.001            0.0007450980392      0.0007450980392")
    
    expect_equal(adjust_indel_rates(rates), output, label="adjusting indel rates")
})
