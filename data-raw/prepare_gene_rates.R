# Gene mutation rates from Nature Genetics 46:944–950 (2014) doi:10.1038/ng.3050

# read-in length of coding sequence of each gene, from Ensembl biomart
url = "http://www.nature.com/ng/journal/v46/n9/extref/ng.3050-S2.xls"

gene_rates = gdata::read.xls(url, sheet="mutation_probabilities")
names(gene_rates)[2] = "hgnc"

gene_rates = gene_rates[, 4:9]  # drop the extra few columns that magically appeared

save(gene_rates, file="data/gene_rates.rda", compress="xz")
