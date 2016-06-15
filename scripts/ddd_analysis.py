"""
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import argparse
import sys

import pandas

from mupit.open_ddd_data import standardise_ddd_de_novos, open_known_genes, \
    get_ddd_rates
from mupit.gene_enrichment import analyse_enrichment
from mupit.write_json_probands import write_probands_by_gene

def get_options():
    """ get the command line options.
    """
    
    parser = argparse.ArgumentParser(description="script to analyse enrichment"
        "of de novo mutations in genes in probands")
    parser.add_argument("--rates",
        help="Path to table of mutation rates.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.mutation_rates.2015-11-24.txt")
    parser.add_argument("--de-novos",
        help="Path to DDD de novo dataset.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-11-24.txt")
    parser.add_argument("--validations",
        help="Path to validation results.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-11-24.txt")
    parser.add_argument("--families",
        help="Path to families PED file.",
        default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt")
    parser.add_argument("--trios",
        help="Path to file listing complete trios.",
        default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/trios.txt")
    parser.add_argument("--known-genes",
        help="path to table of known developmental disorder genes.",
        default="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter")
    parser.add_argument("--diagnosed", help="Path to diagnosed probands file.")
    parser.add_argument("--skip-ddd", default=False, action="store_true",
        help="whether to remove all the DDD probands, and run with the " \
            "external subsets alone.")
    parser.add_argument("--external-cohorts",
        help="Path to table of proband counts in other published de novo datasets.")
    parser.add_argument("--external-variants",
        help="Path to table of de novo mutations in other published de novo datasets.")
    parser.add_argument("--meta-subset",
        help="Comma-separated list of phenotypes eg " \
            "intellectual_disability,autism. This list determines the subset " \
            "of external studies to include. (defaults to using all subsets, " \
            "if the meta-analysis flag is also used).")
    parser.add_argument("--out-manhattan", help="Path to put PDF of manhattan plot.")
    parser.add_argument("--out-probands-by-gene",
        help="Path to put json file of probands per gene.")
    parser.add_argument("--out-enrichment",
        help="Path to put file of enrichment testing results.")
    parser.add_argument("--out-clustering",
        help="Path to put file of enrichment testing results.")
    
    args = parser.parse_args()
    
    if not ((args.external_cohorts is None and args.external_variants is None) | \
            (args.external_cohorts is not None and args.external_variants is not None)):
        sys.exit("You have to either use both of the external-cohorts and "
            "external-variants arguments, or not use either of them.")
    
    return args

def count_trios(families_path, trios_path, diagnosed_path,
        known_genes_path, meta_cohort=None, meta_variants=None, meta_subset=None, skip_ddd=False):
    """ defines the cohort sizes, used to get the overall population size
    
    Args:
        families_path: path to DDD family relationships file, in ped format,
            containing proband IDs and sex information
        trios_path: path to table of probands in complete trios.
        diagnosed_path: path to table of probands with diagnoses
        known_genes_path: path to table of known developmental disorder genes
        meta_cohort: path to table of counts of probands in external exome and
            genome sequencing studies.
        meta_variants: path to table of de novo mutations from external exome
            and genome sequencing studies.
        meta_subset: string of comma-separated list of phenotypes to include in
            the meta-analysis, or None.
        skip_ddd: boolean of whether to not use any DDD probands.
    
    Returns:
        dictionary with total counts of trios with male and female offspring
    """
    
    male = 0
    female = 0
    
    if not skip_ddd:
        ddd_male, ddd_female = count_ddd_trios(families_path, trios_path, diagnosed_path)
        male += ddd_male
        female += ddd_female
    
    if meta_cohort is not None:
        (external_male, external_female) = count_external_trios(meta_cohort,
            meta_variants, known_genes_path, diagnosed_path, meta_subset)
        male += external_male
        female += external_female
    
    return {"male": int(male), "female": int(female)}

def count_ddd_trios(families_path, trios_path, diagnosed_path):
    """ count the male and female probands in the complete DDD trios
    
    Args:
        families_path: path to DDD family relationships file, in ped format,
            containing proband IDs and sex information
        trios_path: path to table of probands in complete trios.
        diagnosed_path: path to table of probands with diagnoses
    
    Returns:
        tuple of male and female proband counts.
    """
    
    # load proband information, then select the proband who have exome sequence
    # available for both parents.
    families = pandas.read_table(families_path, sep="\t")
    trios = pandas.read_table(trios_path, sep="\t")
    proband_ids = trios["proband_stable_id"]
    probands = families[families["individual_id"].isin(proband_ids)]
    
    # get the number of trios studied in our data for each sex
    sex = probands["sex"].value_counts()
    male = sex[["M"]]
    female = sex[["F"]]
    
    if diagnosed_path is not None:
        # remove probands in DDD, unless we are not using the DDD probands.
        diagnosed = pandas.read_table(diagnosed_path, sep="\t")
        diagnosed = diagnosed[~diagnosed[["person_id", "sex"]].duplicated()]
        
        male -= sum(diagnosed["sex"].isin(["Male", "male", "M", "m"]))
        female -= sum(diagnosed["sex"].isin(["Female", "female", "F", "f"]))
    
    return (male, female)

def count_external_trios(meta_cohort, meta_variants, known_genes_path, diagnosed, meta_subset=None):
    """ defines the cohort sizes, used to get the overall population size
    
    Args:
        meta_cohort: path to table of counts of probands in external exome and
            genome sequencing studies.
        meta_variants: path to table of de novo mutations from external exome
            and genome sequencing studies.
        known_genes_path: path to table of known developmental disorder genes
        remove_diagnosed: boolean of whether to remove probands with diagnostic
            variants.
        meta_subset: string of comma-separated list of phenotypes to include in
            the meta-analysis, or None.
    
    Returns:
        tuple of male and female proband counts.
    """
    
    cohorts = pandas.read_table(meta_cohort, sep="\t")
    if meta_subset is not None:
        cohorts = cohorts[cohorts["study_phenotype"].isin(meta_subset.split(","))]
    
    male = sum(cohorts["unique_male"])
    female = sum(cohorts["unique_female"])
    
    if diagnosed is not None:
        variants = pandas.read_table(meta_variants, sep="\t", compression="gzip")
        
        if meta_subset is not None:
            variants = variants[variants["study_phenotype"].isin(meta_subset.split(","))]
        
        known_genes = open_known_genes(known_genes_path)
        diagnosed = variants[variants["hgnc"].isin(known_genes["gene"][known_genes["dominant"]]) |
            ((variants["sex"] == "male") & variants["hgnc"].isin(known_genes["gene"][known_genes["hemizygous"]]))]
        diagnosed = diagnosed[~diagnosed[["person_id", "sex"]].duplicated()]
        
        # decrement for the diagnosed external individuals of each sex
        male -= sum(diagnosed["sex"] == "male")
        female -= sum(diagnosed["sex"] == "female")
    
    return (male, female)

def get_de_novos(de_novos_path, validations_path, diagnosed_path,
    known_genes_path, meta_variants=False, meta_subset=None, skip_ddd=False):
    """ combine datasets listing de novo mutations into a single data frame
    
    Args:
        diagnosed: list of IDs and sex for probands with diagnoses in the DDD
        meta: true/false for whether to include meta-analysis populations
    
    Returns:
        data frame with columns for HGNC, chrom, position, consequence, SNV or
        INDEL type, and study ID.
    """
    
    columns = ["person_id", "sex", "chrom", "start_pos", "end_pos",
        "ref_allele", "alt_allele", "hgnc", "consequence", "study_code",
        "publication_doi", "study_phenotype", "type"]
    variants = pandas.DataFrame(columns=columns)
    
    if not skip_ddd:
        ddd_vars = open_ddd_variants(de_novos_path, validations_path, diagnosed_path)
        variants = variants.append(ddd_vars, ignore_index=True)
    
    if meta_variants:
        external_vars = open_external_variants(meta_variants, meta_subset,
            diagnosed_path, known_genes_path)
        variants = variants.append(external_vars, ignore_index=True)
    
    return variants

def open_ddd_variants(de_novos_path, validations, diagnosed=None):
    variants = standardise_ddd_de_novos(de_novos_path)
    
    # remove diagnosed patients, if maximising power
    if diagnosed is not None:
        diagnosed = pandas.read_table(diagnosed, sep="\t")
        variants = variants[~variants["person_id"].isin(diagnosed["person_id"])]
    
    validations = pandas.read_table(validations, sep="\t")
    variants = variants.merge(validations, how="left",
        on=["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
            "alt_allele", "hgnc", "consequence"])
    
    # drop out the variants that failed to validate (i.e. were false positives,
    # or inherited)
    variants = variants[~variants["status"].isin(["false_positive", "inherited"])]
    del variants["status"]
    
    return variants

def open_external_variants(meta_variants, meta_subset, diagnosed_path, known_genes_path):
    variants = pandas.read_table(meta_variants, sep="\t", compression="gzip")
    if meta_subset is not None:
        variants = variants[variants["study_phenotype"].isin(meta_subset.split(","))]
    
    if diagnosed_path is not None:
        known_genes = open_known_genes(known_genes_path)
        variants = variants[~(variants["hgnc"].isin(known_genes["gene"][known_genes["dominant"]]) |
            ((variants["sex"] == "male") & variants["hgnc"].isin(known_genes["gene"][known_genes["hemizygous"]])))]
    
    return variants

def main():
    args = get_options()
    rates = get_ddd_rates(args.rates)
    
    # analyse the de novos
    trios = count_trios(args.families, args.trios, args.diagnosed,
        args.known_genes, args.external_cohorts, args.external_variants, args.meta_subset, args.skip_ddd)
    de_novos = get_de_novos(args.de_novos, args.validations, args.diagnosed,
        args.known_genes, args.external_variants, args.meta_subset, args.skip_ddd)
    enriched = analyse_enrichment(de_novos, trios,
        plot_path=args.out_manhattan, rates=rates)
    
    # write the enrichment results to a table
    if args.out_enrichment is not None:
        enriched.to_csv(args.out_enrichment, sep="\t", index=False, na_rep="NA")
    
    # and write a list of probands with de novos per gene to a file. This is
    # for HPO similarity testing, so can only be used with DDD samples, since we
    # only have HPO phenotypes available for those individuals.
    if args.external_cohorts is None and args.out_probands_by_gene is not None:
        write_probands_by_gene(de_novos, args.out_probands_by_gene)
    
    if args.out_clustering is not None:
        # write the set of de novos for clustering analysis
        de_novos[["hgnc", "chrom", "start_pos", "consequence", "type"]].to_csv(
            args.out_clustering, sep="\t", index=False)

if __name__ == '__main__':
    main()
