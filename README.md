### Mupit: de novo mutation recurrence significance testing
program to calculate the significance of seeing N DNMs of a specific
combination of functional types in a particular gene in M trios

RATIONALE: use gene coding sequence to predict rate of DNMs in coding sequence
for each gene of different functional classes, then estimate the probability
of seeing the observed combination of different functional classes of DNMs
assuming number of DNMs in each class is Poisson distributed

initial implementation: use genome-wide mutation rate and scale by length of
coding sequence, use genome-wide average of functional consequences of coding
mutations from Kryukov et al 2007


#### POTENTIAL FUTURE IMPROVEMENTS (highest priority first)

- [x] adapt for chrX
- [x] FDR estimates
- [x] analytical rather than permutation
- [ ] use actual number of exons to predict essential_splice_site mutations
- [x] predict mutation rates from base composition
- [ ] use estimate of de novo mutation discovery power in a gene to 
    better estimate gene-specific mutation rate.
- [x] get CDS length for all genes. use longest transcript if >1
- [ ] calculate coding sequence length according to intersection of exome 
    targeted regions and the union of all transcripts for a gene.
- [ ] account for incomplete sensitivity for DNMs, especially indels
- [x] look at clustering of de novos in genes with recurrent mutations 

#### Format
Read in table of genes and the numbers of families with SNV DNMs
within different functiona classes and indel DNMS in different functional
classes, output same table with added column of probability of seeing that
combination of DNMs.

#### Input data
validated de novo mutations in TSV format:

proband ID | chrom | start | stop | ref | alt | HGNC symbol | VEP consequence | study | DOI | phenotype | SNV or indel

note: current test input file is numbers of mutations, not number of families,
some families have >1 mutation in the same gene
