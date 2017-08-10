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

from mupit.open_ddd_data import open_known_genes, standardise_ddd_de_novos

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--ddd-1k-diagnoses", help="Path to clinically reviewed diagnoses in the DDD 1K.",
        default="/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx")
    parser.add_argument("--updated-diagnoses",
        help="Path to table of clinially reviewed diagnoses in the remainder.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novo_data/ddd.reported_variants.txt")
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
        
    return parser.parse_args()

def get_ddd_diagnostic_cnvs(path):
    """ load the CNVs
    """
    
    reviewed_cnvs = pandas.read_excel(path, sheetname="CNVs reviewed")
    
    # redefine the columns we want from the CNVs table
    cnvs = reviewed_cnvs[reviewed_cnvs["Diagnostic?"] > 0].copy()
    cnvs["person_id"] = cnvs["DDDP_ID"]
    cnvs["chrom"] = cnvs["Chr"].astype(str)
    cnvs["start_pos"] = cnvs["Start"].astype('object')
    cnvs["end_pos"] = cnvs["Stop"].astype('object')
    cnvs["ref_allele"] = None
    cnvs["alt_allele"] = None
    cnvs["hgnc"] = None
    cnvs["inheritance"] = cnvs["Inheritance"]
    cnvs['confirmed'] = 'yes'
    cnvs["type"] = "cnv"
    
    return cnvs[["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "confirmed", "type"]]

def get_ddd_diagnostic_snvs(path):
    """ load the diagnostic SNVs
    """
    
    reviewed_snvs = pandas.read_excel(path, sheetname="Exome variants reviewed")
    snvs = reviewed_snvs[reviewed_snvs["DECISION"] == "Yes"].copy()
    
    snvs["person_id"] = snvs["proband"]
    snvs["start_pos"] = snvs["position"].astype('object')
    snvs['chrom'] = snvs['chrom'].astype(str)
    
    alleles = snvs["ref/alt_alleles"].str.split("/")
    snvs["ref_allele"] = alleles.str.get(0)
    snvs["alt_allele"] = alleles.str.get(1)
    snvs["end_pos"] = snvs["start_pos"] + snvs["ref_allele"].str.len() - 1
    snvs["end_pos"] = snvs["end_pos"].astype('object')
    
    snvs["hgnc"] = snvs["gene"]
    snvs.loc[:, "inheritance"] = snvs["inheritance (DECIPHER compatible)"]
    snvs['confirmed'] = 'yes'
    
    # determine the type of variant
    is_indel = (snvs['ref_allele'].str.len() > 1) | (snvs['alt_allele'].str.len() > 1)
    snvs['type'] = is_indel.map({True: 'indel', False: 'snv'})
    
    return snvs[["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "confirmed", "type"]]

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
    other['confirmed'] = 'yes'
    
    other['chrom'] = other['chrom'].astype(str)
    other['start_pos'] = other['start_pos'].astype('object')
    other['end_pos'] = other['end_pos'].astype('object')
    
    recode = {"Mosaicism": "mosaic_cnv", "UPD": "uniparental_disomy"}
    other["type"] = other["type"].map(recode)
    
    return other[["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "confirmed", "type"]]

def get_updated_diagnoses(path):
    ''' load clinically reviewed diagnoses made in later DDD datafreezes
    '''
    
    data = pandas.read_table(path)
    
    data['start_pos'] = data['start'].astype('object')
    data['end_pos'] = data['end'].astype('object')
    data['hgnc'] = data['symbol']
    data['chrom'] = data['chrom'].astype(str)
    
    # subset to the pathogenic variants
    required = ['Definitely pathogenic', 'Likely pathogenic']
    data = data[data['pathogenicity'].isin(required)]
    data['confirmed'] = 'yes'
    
    return data[["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
        "alt_allele", "hgnc", "inheritance", "confirmed", "type"]]

def get_reviewed(path, families_path, updated_path):
    """ find clinically reviewed diagnoses, to exclude probands from analyses
    
    Args:
        path: path to file defining diagnosed probands
    
    Returns:
        A list containing vectors with DDD IDs, and sex of the diagnosed
        probands.
    """
    
    cnvs = get_ddd_diagnostic_cnvs(path)
    snvs = get_ddd_diagnostic_snvs(path)
    other = get_other_diagnostic_variants(path)
    updated = get_updated_diagnoses(updated_path)
    
    diagnosed = pandas.concat([cnvs, snvs, other, updated], axis=0, ignore_index=True)
    
    # remove duplicate diagnoses
    is_dup = diagnosed.duplicated(subset=['person_id', 'hgnc'], keep='first')
    diagnosed = diagnosed[~is_dup]
    
    # relabel the inheritance types
    recode = {'DNM': 'de_novo',
        'De novo constitutive': 'de_novo',
        'De novo mosaic': 'de_novo',
        'Biparental': 'biparental',
        'Maternal': 'maternal',
        '[Pp]aternally inherited, constitutive in father': 'paternal',
        '[Mm]aternally inherited, constitutive in mother ': 'maternal',
        'Maternally inherited, mosaic in mother': 'de_novo',
        'Paternally inherited, mosaic in father': 'de_novo',}
    diagnosed['inheritance'] = diagnosed['inheritance'].map(recode)
    
    # get the sex info for each proband
    families = pandas.read_table(families_path, sep="\t")
    families['sex'] = families['sex'].map({'M': 'male', 'F': 'female'})
    recode = dict(zip(families['individual_id'], families['sex']))
    
    diagnosed['sex'] = diagnosed['person_id'].map(recode)
    
    return diagnosed[["person_id", "sex", "chrom", "start_pos", "end_pos",
         "ref_allele", "alt_allele", "hgnc", "inheritance", "confirmed", "type"]]

def get_low_pp_dnm_validations(path):
    """ load the validation data for candidates with low pp_dnm scores
    """
    
    data = pandas.read_excel(path, sheetname="Sheet1")
    
    data.rename(columns={'pos': 'start_pos'}, inplace=True)
    
    data["chrom"] = data["chrom"].astype(str)
    is_indel = (data['ref_allele'].str.len() > 1) | (data['alt_allele'].str.len() > 1)
    data['type'] = is_indel.map({True: 'indel', False: 'snv'})
    
    return data[["person_id", "chrom", "start_pos", "ref_allele",
        "alt_allele", "type", "status"]]

def get_current_de_novos(path):
    """ load the current set of de novos to be analysed
    """
    
    variants = standardise_ddd_de_novos(path, extra_columns=['pp_dnm'])
    variants["inheritance"] = "de_novo"
    variants['confirmed'] = "no"
    
    return variants[["person_id", "sex", "chrom", "start_pos", "end_pos",
         "ref_allele", "alt_allele", "hgnc", "inheritance", "confirmed", "type",
         "consequence", "pp_dnm"]]

def check_for_match(site, initial):
    """ checks if a sites has a match in a previous dataset
    
    Some de novos in the current dataset were also present in a previous
    datafreeze. We need to remove these so as to not double count diagnostic
    variants. Unfortunately, some sites have shifted location slightly such as
    indels, which can be difficult to locate.
    """
    
    data = initial[(initial["person_id"] == site["person_id"])]
    data = data[(data["chrom"] == site["chrom"]) & ~data["chrom"].isnull()]
    
    if len(data) == 0:
        return False
    
    delta = abs(data["start_pos"] - site["start_pos"])
    
    if min(delta) > 20:
        return False
    
    return sum(delta == min(delta)) == 1

def get_diagnosed(diagnosed_path, updated_path, de_novo_path,
        low_pp_dnm_validations_path, known_genes_path, families_path,
        recessive_path):
    """ find probands likely to have diagnoses, to exclude them from our data
    
    Args:
        path: path to file defining diagnosed probands
    
    Returns:
        A table of probands with diagnoses
    """
    
    initial_diagnosed = get_reviewed(diagnosed_path, families_path, updated_path)
    
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
    
    # define the sufficiently pathogenic consequences
    permitted = ["missense_variant", "frameshift_variant", "stop_gained",
        "splice_donor_variant", "splice_acceptor_variant", "inframe_deletion",
         "conserved_exon_terminus_variant", "initiator_codon_variant",
         "inframe_insertion"]
    
    # remove the nonfunctional variants
    likely_diagnostic = likely_diagnostic[likely_diagnostic["consequence"].isin(permitted)]
    likely_diagnostic = likely_diagnostic[["person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "inheritance", "confirmed", "type"]]
    
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
    
    diagnosed = get_diagnosed(args.ddd_1k_diagnoses, args.updated_diagnoses,
        args.de_novos, args.low_pp_dnm, args.known_genes, args.families,
        args.recessive_diagnoses)
    
    diagnosed.to_csv(args.out, sep="\t", na_rep="NA", index=False)

if __name__ == '__main__':
    main()
