# construct R package data frame of sex-specific numbers of affected probands
# reported on in published exome sequencing of children with developmental 
# disorders (autism, congenital heart disease, intellectual disability, 
# epilepsy, and schizophrenia). 
# 
# This provides the total population of trio-based probands investigated within
# these studies, that is, the population for whom a de novo mutation might have
# been identified. This determines the number of mutations we expect to find
# within the total population.
# 
# The table of sex-counts has been carefully extracted from the publications, 
# accounting for overlap in samples between different publications.
cohorts = read.table(file.path("data-raw", "cohorts.txt"), sep="\t", 
    header=TRUE, stringsAsFactors=FALSE)

save(cohorts, file="data/cohorts.rda", compress="xz")

