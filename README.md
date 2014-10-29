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
- adapt for chrX, lower mutation rate, and number of transmissions required.
      Assume male/female of proband = 1:1, unless known [DONE]
- incorporate FDR estimation [DONE]
- convert to analytical approach from permutation approach [DONE]
- use actual number of exons to predict essential_splice_site mutations,
- use base composition of coding sequence to predict gene-specific mutation
      rate for each class of mutation, rather than just scaling the
      genome-average.
- use estimate of de novo mutation discovery power in a gene to better
      estimate gene-specific mutation rate
- get CDS length for all genes, not just this subset, 455/478 in test data.
      use longest transcript if >1 [DONE]
- calculate coding sequence length according to intersection of exome
      targeted regions and the union of all transcripts for a gene.
- account for incomplete sensitivity for DNMs, especially indels
- look at clustering of de novos in genes with recurrent mutations, within
      protein space

format: read in table of genes and the numbers of families with SNV DNMs
within different functiona classes and indel DNMS in different functional
classes, output same table with added column of probability of seeing that
combination of DNMs.

input data (validated SNV DNMs in TSV format, HGNC_ID, NUM_LOF, NUM_MISSENSE,
NUM_LOF_INDEL, NUM_MISSENSE_INDEL):

note: current test input file is numbers of mutations, not number of families,
some families have >1 mutation in the same gene
