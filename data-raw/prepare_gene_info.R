# Ensembl derived chromosome and coding sequence lengths for all HGNC symbols

# read-in length of coding sequence of each gene, from Ensembl biomart
gene_info = read.delim(file.path("data-raw", "CDS_LENGTH_B37_chr.txt"), header=TRUE)
gene_info$cds_length = gene_info$CDS_LENGTH
gene_info$hgnc = gene_info$ID
gene_info$chrom = gene_info$chr

# only include the renamed columns
gene_info = subset(gene_info, select = c(hgnc, chrom, cds_length))

save(gene_info, file="data/gene_info.rda")

