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
from datetime import datetime

import pandas

from mupit.open_ddd_data import open_known_genes

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--ddd-1k-diagnoses", help="Path to DDD 1K diagnoses.",
        default="/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx")
    parser.add_argument("--de-novos", help="Path to DDD de novo dataset.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-10-12.txt")
    parser.add_argument("--low-pp-dnm", help="Path to low PP_DNM validations.",
        default="/nfs/users/nfs_j/jm33/de_novos.ddd_4k.validation_results.low_pp_dnm.2015-10-02.xlsx")
    parser.add_argument("--recessive-diagnoses", help="Path to additional recessive diagnoses.")
    parser.add_argument("--families", help="Path to families PED file.",
        default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt")
    parser.add_argument("--known-genes", help="Path to table of known developmental disorder genes.",
        default="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter")
    parser.add_argument("--out", help="Path to output file.",
        default="/lustre/scratch113/projects/ddd/users/jm33/ddd_4k.diagnosed.{}.txt".format(datetime.today().strftime("%Y-%m-%d")))
        
    args = parser.parse_args()
    
    return args

def get_ddd_diagnostic_cnvs(path):
    """ load the CNVs
    """
    
    reviewed_cnvs = pandas.read_excel(path, sheetname="CNVs reviewed")
    
    # redefine the columns we want from the CNVs table
    cnvs = reviewed_cnvs[reviewed_cnvs["Diagnostic?"] > 0].copy()
    cnvs["person_id"] = cnvs["DDDP_ID"]
    cnvs["chrom"] = cnvs["Chr"]
    cnvs["start_pos"] = cnvs["Start"]
    cnvs["end_pos"] = cnvs["Stop"]
    cnvs["ref_allele"] = None
    cnvs["alt_allele"] = None
    cnvs["hgnc"] = None
    cnvs["inheritance"] = cnvs["Inheritance"]
    cnvs["type"] = "cnv"
    
    cnvs = cnvs[["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "type"]]
    
    return cnvs

def get_ddd_diagnostic_snvs(path):
    """ load the diagnostic SNVs
    """
    
    reviewed_snvs = pandas.read_excel(path, sheetname="Exome variants reviewed")
    snvs = reviewed_snvs[reviewed_snvs["DECISION"] == "Yes"].copy()
    
    snvs["person_id"] = snvs["proband"]
    snvs["start_pos"] = snvs["position"]
    
    alleles = snvs["ref/alt_alleles"].str.split("/")
    snvs["ref_allele"] = alleles.str.get(1)
    snvs["alt_allele"] = alleles.str.get(2)
    snvs["end_pos"] = snvs["start_pos"] + snvs["ref_allele"].str.len() - 1
    
    snvs["hgnc"] = snvs["gene"]
    snvs.loc[:, "inheritance"] = snvs["inheritance (DECIPHER compatible)"]
    
    # determine the type of variant
    snvs["type"] = "snv"
    snvs["type"][(snvs["ref_allele"].str.len() > 1) | (snvs["alt_allele"].str.len() > 1)] = "indel"
    
    snvs = snvs[["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "type"]].copy()
    
    return snvs

def get_other_diagnostic_variants(path):
    """ now load the other diagnostic variants
    """
    
    other = pandas.read_excel(path, sheetname="UPD&Mosaicism", header=None)
    
    other.columns = ["person_id", "type"]
    other["chrom"] = None
    other["start_pos"] = None
    other["end_pos"] = None
    other["ref_allele"] = None
    other["alt_allele"] = None
    other["hgnc"] = None
    other["inheritance"] = None
    
    other["type"] = other["type"].str.replace("Mosaicism", "mosaic_cnv")
    other["type"] = other["type"].str.replace("UPD", "uniparental_disomy")
    
    other = other[["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "type"]]
    
    return other

def get_previous(path, families_path):
    """ find diagnosed probands in the DDD study, to exclude them from our data
    
    Args:
        path: path to file defining diagnosed probands
    
    Returns:
        A list containing vectors with DDD IDs, and sex of the diagnosed
        probands.
    """
    
    cnvs = get_ddd_diagnostic_cnvs(path)
    snvs = get_ddd_diagnostic_snvs(path)
    other = get_other_diagnostic_variants(path)
    
    diagnosed = pandas.concat([cnvs, snvs, other], axis=0, ignore_index=True)
    diagnosed["chrom"] = diagnosed["chrom"].astype(str)
    
    # relabel the inheritance types
    inheritance = diagnosed["inheritance"]
    inheritance[inheritance.str.match("DNM") & ~inheritance.isnull()] = "de_novo"
    inheritance[inheritance.str.match("De novo constitutive") & ~inheritance.isnull()] = "de_novo"
    inheritance[inheritance.str.match("Maternal") & ~inheritance.isnull()] = "maternal"
    inheritance[inheritance.str.match("Biparental") & ~inheritance.isnull()] = "biparental"
    inheritance[inheritance.str.match("[Pp]aternally inherited, constitutive in father") & ~inheritance.isnull()] = "paternal"
    inheritance[inheritance.str.match("[Mm]aternally inherited, constitutive in mother") & ~inheritance.isnull()] = "maternal"
    diagnosed["inheritance"] = inheritance
    
    # get the sex info for each proband
    families = pandas.read_table(families_path, sep="\t")
    families = families[["individual_id", "sex"]].copy()
    families.columns = ["person_id", "sex"]
    
    families["sex"] = families["sex"].str.replace("M", "male")
    families["sex"] = families["sex"].str.replace("F", "female")
    
    diagnosed = diagnosed.merge(families, how="left", on="person_id")
    
    return diagnosed

def get_low_pp_dnm_validations(path):
    """ load the validation data for candidates with low pp_dnm scores
    """
    
    low_pp_dnm = pandas.read_excel(path, sheetname="Summary_Final_forDB")
    
    # only select the most useful columns
    low_pp_dnm = low_pp_dnm[["ID", "CHR", "POS", "REF", "ALT", "TYPE",
        "manual_score"]]
        
    low_pp_dnm.columns = ["person_id", "chrom", "start_pos", "ref_allele",
        "alt_allele", "type", "status"]
    
    low_pp_dnm["chrom"] = low_pp_dnm["chrom"].astype(str)
    
    # recode the validation status codes to more understandable codes
    status = low_pp_dnm["status"]
    status[status == "dnm"] = "de_novo"
    status[status == "dnm_low_alt"] = "de_novo"
    status[status == "fp"] = "false_positive"
    status[status == "inherited_pat"] = "inherited"
    status[status == "p/u"] = "uncertain"
    status[status == "unclear"] = "uncertain"
    status[status == "parental_mosaic"] = "de_novo"
    low_pp_dnm["status"] = status
    
    low_pp_dnm["type"][low_pp_dnm["type"] == "DENOVO-SNP"] = "snv"
    low_pp_dnm["type"][low_pp_dnm["type"] == "DENOVO-INDEL"] = "indel"
    
    return low_pp_dnm

def get_current_de_novos(path):
    """ load the current set of de novos to be analysed
    """
    
    variants = pandas.read_table(path, sep="\t")
    # standardise the SNV or INDEL flag
    variants["type"] = "indel"
    variants["type"][variants["var_type"] == "DENOVO-SNP"] = "snv"
    
    # standardise the columns, and column names
    variants["person_id"] = variants["person_stable_id"]
    variants["chrom"] = variants["chrom"].astype(str)
    variants["start_pos"] = variants["pos"]
    variants["ref_allele"] = variants["ref"]
    variants["alt_allele"] = variants["alt"]
    variants["end_pos"] = variants["start_pos"] + variants["ref_allele"].str.len() - 1
    
    variants["hgnc"] = variants["symbol"]
    variants["inheritance"] = "de_novo"
    
    variants["sex"] = variants["sex"].str.replace("M", "male")
    variants["sex"] = variants["sex"].str.replace("F", "female")
    
    return variants

def check_for_match(site, initial, pos=False):
    """ checks if a sites has a match in a previous dataset
    
    Some de novos in the current dataset were also present in a previous
    datafreeze. We need to remove these so as to not double count diagnostic
    variants. Unfortunately, som sites have shifted location slightly such as
    indels, which can be difficult to locate.
    """
    
    # set the missing match type dependent on whether we are looking for whether
    # we have a match, or the matching position.
    if pos:
        no_match = None
    else:
        no_match = False
    
    vars = initial[(initial["person_id"] == site["person_id"])]
    vars = vars[(vars["chrom"] == site["chrom"]) & ~vars["chrom"].isnull()]
    
    if len(vars) == 0:
        return no_match
    
    delta = abs(vars["start_pos"] - site["start_pos"])
    
    if min(delta) > 20:
        return no_match
    
    close = delta == min(delta)
    
    if sum(close) == 1:
        if not pos:
            return True
        else:
            return vars["start_pos"][close]
    else:
        return no_match

def get_ddd_diagnosed(diagnosed_path, de_novo_path, low_pp_dnm_validations_path,
        known_genes_path, families_path, recessive_path):
    """ find probands likely to have diagnoses, to exclude them from our data
    
    Args:
        path: path to file defining diagnosed probands
    
    Returns:
        A table of probands with diagnoses
    """
    
    initial_diagnosed = get_previous(diagnosed_path, families_path)
    
    known_genes = open_known_genes(known_genes_path)
    variants = get_current_de_novos(de_novo_path)
    
    # the candidates with low pp_dnm (< 0.9) in DDG2P genes were attempted to
    # validate. Those that did, we can swap the pp_dnm to 1, since these are now
    # high confidence de novos
    low_pp_dnm = get_low_pp_dnm_validations(low_pp_dnm_validations_path)
    variants = variants.merge(low_pp_dnm[["person_id", "chrom", "start_pos", "status"]],
        how="left", on=["person_id", "chrom", "start_pos"])
    variants["pp_dnm"][variants["status"] == "de_novo"] = 1
    
    # get the set of de novos from the current dataset that are likely to be
    # diagnostic. These are de novos in genes with dominant modes of inheritance,
    # or chrX de novos in males in genes with hemizygous mode of inheritance,
    # and where the site is high confidence (as determined by having a high
    # pp_dnm, or a missing pp_dnm)
    dominant = known_genes["gene"][known_genes["dominant"]]
    hemizygous = known_genes["gene"][known_genes["hemizygous"]]
    likely_diagnostic = variants[(variants["hgnc"].isin(dominant) | \
        (variants["hgnc"].isin(hemizygous) & (variants["sex"] == "male")))
        & ((variants["pp_dnm"] > 0.1) | variants["pp_dnm"].isnull())]
    
    # remove the synonymous variants
    likely_diagnostic = likely_diagnostic[likely_diagnostic["consequence"] != "synonymous_variant"]
    likely_diagnostic = likely_diagnostic[["person_id", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "inheritance", "type",
        "sex"]]
    
    # remove the sites from the likely diagnoses that form part of a confirmed
    # diagnosis
    in_prev = likely_diagnostic.apply(axis=1, func=check_for_match, initial=initial_diagnosed)
    likely_diagnostic = likely_diagnostic[~in_prev]
    
    diagnosed = initial_diagnosed.append(likely_diagnostic, ignore_index=True)
    
    if recessive_path is not None:
        recessive = pandas.read_table(recessive_path, sep="\t")
        diagnosed = diagnosed.append(recessive, ignore_index=True)
    
    return diagnosed

def main():
    
    args = get_options()
    
    diagnosed = get_ddd_diagnosed(args.ddd_1k_diagnoses, args.de_novos,
        args.low_pp_dnm, args.known_genes, args.families, args.recessive_diagnoses)
    
    diagnosed.to_csv(args.out, sep="\t", na_rep="NA", index=False)

if __name__ == '__main__':
    main()
