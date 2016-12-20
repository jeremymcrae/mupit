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

import math

import numpy
import pandas
from scipy.stats import chi2

def fishersMethod(x):
    """ function to combine p values, using Fisher's method
    
    Args:
        x: list of P-values for a gene
    
    Returns:
        combined P-value
    """
    
    x = [ val for val in x if not math.isnan(val) ]
    
    if len(x) == 0:
        return numpy.nan
    
    return chi2.sf(-2 * sum(numpy.log(x)), 2 * len(x))

def combine_enrichment_and_clustering(enriched, clust):
    """ combine P values from enrichment and clustering tests into a single P value
     
    Args:
        enriched: dataframe of de novo enrichment results
        clust: dataframe of de novo clustering results
     
    Returns:
        a merged dataset where the P values have been combined
    """
    
    # read in p values from clustering analysis, only for genes with >1 mutation
    clust = pandas.pivot_table(clust, index=["gene_id"],
        columns=["mutation_category"], values="probability", aggfunc=numpy.mean)
    
    clust["hgnc"] = list(clust.index)
    columns = ["missense", "nonsense"]
    rename = dict(zip(columns, [ "p_{}_clust".format(x) for x in columns ]))
    clust = clust.rename(columns=rename)
    
    # merge the datasets
    merged = enriched.merge(clust, how="left", on="hgnc")
    
    # calculate a combined p-value for each gene. We don't expect the
    # loss-of-function de novos to be clustered, so we don't use that.
    p_values = merged[["p_func", "p_missense_clust"]]
    merged["p_combined"] = p_values.apply(fishersMethod, axis=1)
    
    # calculate minimum p value across LoF and func + clustering tests
    merged["p_min"] = merged[["p_lof", "p_combined"]].min(axis=1)
    
    return merged

def combine_tests(meta_clust, meta_enrich, clust, enrich, pheno_path=None):
    """ find the most significant P value for each gene from the P values from
    different subsets and different tests.
    
    Args:
        meta_clust: path to clustering results for the meta-analysis subset
        meta_enrich: path to enrichment results for the meta-analysis subset
        clust: path to clustering results for the ddd only subset
        enrich: path to enrichment results for the ddd only subset
        pheno_path: path to phenotype similarity testing results
    
    Returns:
        data frame with the columns from all the datasets, as well as minimum
        P values from each subset for each gene, and overall minimum P values
        for each gene.
    """
    
    # load all the data files
    clust = pandas.read_table(clust, sep="\t")
    enrich = pandas.read_table(enrich, sep="\t")
    meta_clust = pandas.read_table(meta_clust, sep="\t")
    meta_enrich = pandas.read_table(meta_enrich, sep="\t")
    
    meta = combine_enrichment_and_clustering(meta_enrich, meta_clust)
    ddd = combine_enrichment_and_clustering(enrich, clust)
    
    # if we have phenotypic similarity results, merge them with the other results
    if pheno_path is not None:
        phenotypes = pandas.read_table(pheno_path, sep="\t")
        ddd = ddd.merge(phenotypes, how="outer", on="hgnc")
    
    # need names that are more informative as same across files, add prefix
    columns = ["lof_indel", "lof_snv", "missense_indel", "missense_snv",
        "p_lof", "p_func", "p_missense_clust", "p_nonsense_clust", "gene_id",
        "p_combined", "hpo_similarity_p_value", "p_min"]
    ddd = ddd.rename(columns=dict(zip(columns, [ "ddd.{}".format(x) for x in columns ])))
    meta = meta.rename(columns=dict(zip(columns, [ "meta.{}".format(x) for x in columns ])))
    
    # merge together files, focusing on genes with DNMs in DDD
    merged = meta.merge(ddd, how="outer", on=["hgnc", "chrom"])
    merged["p_min"] = merged[["ddd.p_min", "meta.p_min"]].min(axis=1)
    
    return merged
