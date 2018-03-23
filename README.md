[![Build Status](https://travis-ci.org/jeremymcrae/mupit.svg?branch=master)](https://travis-ci.org/jeremymcrae/mupit)

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

#### Install
```sh
pip install git+git://github.com/jeremymcrae/mupit.git --user

# Alternatively:
git clone https://github.com/jeremymcrae/mupit.git
cd mupit
python setup.py install --user
```

#### Usage (in python):
```python
from mupit.gene_enrichment import analyse_enrichment
import pandas

trios = {"female": 0, "male": 1}
de_novos = pandas.DataFrame({
    "person_id": ["temp"],
    "chrom": ["16"],
    "start_pos": [89348744],
    "end_pos": [89348744],
    "ref_allele": ["A"],
    "alt_allele": ["G"],
    "hgnc": ["ANKRD11"],
    "consequence": ["missense_variant"],
    "study_code": [None],
    "publication_doi": [None],
    "study_phenotype": [None],
    "type": ["snv"],
    })

enrichment = analyse_enrichment(de_novos, trios)
```

You can also pass in in a user-specified set of log10-scaled gene-based mutation
rates with a dataframe such as:

 hgnc  | length |  mis  |  non  |  splice_site  |  syn  | frameshift
-------|--------|-------|-------|---------------|-------|-----------
 AADAC |  4364  | -4.98 | -6.26 | -6.74         | -5.37 | -6.08
ARID1B |  6869  | -4.06 | -5.42 | -6.06         | -4.30 | -5.33
 KMT2A | 11918  | -3.95 | -5.12 | -5.80         | -4.34 | -5.09
 SCN2A |  6109  | -4.26 | -5.47 | -5.30         | -4.69 | -5.38
 
The columns are HGNC symbol, CDS length, rate of missense mutations, rate of
nonsense mutations, rate of canonical splice site mutations, rate of synonymous
mutations, rate of frameshift mutations. The header needs to use the column
names given in the table above (hgnc, length, mis, non, css, syn, frameshift).

The gene rates can be used with:
```python
rates = pandas.DataFrame({
    "hgnc": ["ANKRD11"],
    "chrom": ["16"],
    "mis": [9e-5],
    "non": [1e-6],
    "splice_site": [1e-7],
    "syn": [1e-5],
    "frameshift": [1e-6],
    })

enrichment = analyse_enrichment(de_novos, trios, rates=rates)
```

#### Format
Read in table of genes and the numbers of families with SNV DNMs
within different functional classes and indel DNMS in different functional
classes, output same table with added column of probability of seeing that
combination of DNMs.

#### Input data
validated de novo mutations in TSV format:

proband ID | sex | chrom | start | stop | ref | alt | HGNC symbol | VEP consequence | study | DOI | phenotype | SNV or indel

note: current test input file is numbers of mutations, not number of families,
some families have >1 mutation in the same gene
