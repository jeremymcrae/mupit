[![Build Status](https://travis-ci.org/jeremymcrae/mupit.svg?branch=master)]
(https://travis-ci.org/jeremymcrae/mupit)

### Mupit: de novo mutation recurrence significance testing
program to calculate the significance of seeing N DNMs of a specific
combination of functional types in a particular gene in M trios

RATIONALE: use gene coding sequence to predict rate of DNMs in coding sequence
for each gene of different functional classes, then estimate the probability
of seeing the observed combination of different functional classes of DNMs
assuming number of DNMs in each class is Poisson distributed

initial implementation: use genome-wide mutation rate and scale by length of
coding sequence, use genome-wide average of functional consequences of coding
mutations from [Kryukov et al. 2007](http://dx.doi.org/10.1086%2F513473).

#### Usage (requires >= R 3.1.0) with:
```R
# obtain and install
library(devtools) # if necessary install with install.packages("devtools")
devtools::install_github("jeremymcrae/mupit")

# load the package
library(mupit)

# get the counts of male and female probands
trios = list(males=1, females=0)
de_novos = data.frame(person_id="temp", chrom=16, start_pos=89348744, end_pos=89348744,
    ref_allele="A", alt_allele="G", hgnc="ANKRD11", consequence="missense_variant",
    study_code=NA, publication_doi=NA, study_phenotype=NA, type="snv")
enriched = analyse_gene_enrichment(de_novos, trios)
```

You can also pass in in a user-specified set of log10-scaled gene-based mutation
rates with a dataframe such as:

 hgnc  | length |  mis  |  non  |  css  |  syn  | frameshift
-------|--------|-------|-------|-------|-------|-----------
 AADAC |  4364  | -4.98 | -6.26 | -6.74 | -5.37 | -6.08
ARID1B |  6869  | -4.06 | -5.42 | -6.06 | -4.30 | -5.33
 KMT2A | 11918  | -3.95 | -5.12 | -5.80 | -4.34 | -5.09
 SCN2A |  6109  | -4.26 | -5.47 | -5.30 | -4.69 | -5.38
 
The columns are HGNC symbol, CDS length, rate of missense mutations, rate of
nonsense mutations, rate of canonical splice site mutations, rate of synonymous
mutations, rate of frameshift mutations. The header needs to use the column
names given in the table above (hgnc, length, mis, non, css, syn, frameshift).

The gene rates can thereafter be used with:
```R
enriched = analyse_gene_enrichment(de_novos, trios, rates=RATES_DATAFRAME)
```

#### Potential future improvements (highest priority first)
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

proband ID | sex | chrom | start | stop | ref | alt | HGNC symbol | VEP consequence | study | DOI | phenotype | SNV or indel

note: current test input file is numbers of mutations, not number of families,
some families have >1 mutation in the same gene
