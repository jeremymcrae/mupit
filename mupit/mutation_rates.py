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

from __future__ import division

import tempfile
import urllib

import pandas

from mupit.gtf import convert_gtf
from mupit.util import is_url

def get_default_rates(rates_url="http://www.nature.com/ng/journal/v46/n9/extref/ng.3050-S2.xls",
        gencode_url="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"):
    """ obtain the table of mutation rates from the Samocha et al paper
    
    Rates for all genes can be obtained from the supplementary material of
    Samocha et al. Nature Genetics 2014 doi:10.1038/ng.3050
    
    Args:
        rates_url: url to supplementary mutation rates table
        gencode_url: url to gencode, or local path. This is required to identify
            chromosomes for the genes in the rates data, since we need to know
            the chromosome in order to corrrect rates on chromosome X.
    
    Returns:
        dataframe of mutation rates, with an extra column for summed lof rate
    """
    
    rates = pandas.read_excel(rates_url, sheetname="mutation_probabilities")
    
    # convert rates from log-scaled values, so we can later multiply by the
    # number of transmissions
    columns = ["syn", "mis", "non", "splice_site", "frameshift"]
    rates[columns] = 10 ** rates[columns]
    
    # sort out the required columns and names.
    rates["hgnc"] = rates["gene"]
    gencode = load_gencode(gencode_url)
    recode = dict(zip(gencode["hgnc"], gencode["chrom"]))
    rates["chrom"] = rates["hgnc"].map(recode)
    
    rates = rates[["hgnc", "chrom", "syn", "mis", "splice_site", "frameshift", "non"]]
    
    return rates

def load_gencode(path):
    """ load gencode table with HGNC symbols and chromosome coordinates
    
    Args:
        path: path to gzipped tab-separated table of gencode gene entries. This
            can be either a url, or local path.
    
    Returns:
        pandas dataframe of HGNC symbols and genome coordiantes
    """
    
    gencode = convert_gtf(path)
    
    # restrict outselves to protein coding genes (or genes which are protein
    # coding in at least some individuals)
    gencode = gencode[gencode["gene_type"].isin(["protein_coding",
        "polymorphic_pseudogene"])]
    
    gencode = gencode[gencode["feature"] == "gene"]
    
    # get the required column names, and strip out all unnecessary columns
    gencode["hgnc"] = gencode["gene_name"]
    gencode["chrom"] = [ x.strip("chr") for x in gencode["seqname"].astype(str) ]
    
    gencode = gencode[["hgnc", "chrom", "start", "end"]].copy()
    
    return gencode

def get_expected_mutations(rates, male, female):
    """ gets numbers of expected mutation per gene
    
    Loads gene-based mutation rates, in order to determine the expected number
    of mutations per gene, given the number of studied probands and adjusts for
    sex-chromosome transmissions.
    
    This defaults to the gene-based mutation rates from Nature Genetics
    46:944-950 (2014) doi:10.1038/ng.3050, but we can pass in other gene-based
    mutation rate datasets.
    
    Args:
        rates: pandas dataframe containing per-gene mutation
        male: number of male probands in the dataset
        female: number of female probands in the dataset
    
    Returns:
        a dataframe of mutation rates for genes under different mutation
        classes.
    """
    
    if rates is None:
        rates = get_default_rates()
    
    autosomal = 2 * (male + female)
    
    expected = rates[["hgnc", "chrom"]].copy()
    
    # get the number of expected mutations, given the number of transmissions
    expected["lof_indel"] = rates["frameshift"] * autosomal
    expected["lof_snv"] = (rates[["non", "splice_site"]].sum(axis=1, skipna=True)) * autosomal
    expected["missense_indel"] = (rates["frameshift"] / 9) * autosomal
    expected["missense_snv"] = rates["mis"] * autosomal
    expected["synonymous_snv"] = rates["syn"] * autosomal
    
    # correct for the known ratio of indels to nonsense, and for transmissions
    # on the X-chromosome
    expected = adjust_indel_rates(expected)
    expected = correct_for_x_chrom(expected, male, female)
    
    # subset to the columns we need to estimate enrichment probabilities
    expected = expected[["hgnc", "chrom", "lof_indel", "lof_snv",
        "missense_indel", "missense_snv", "synonymous_snv"]]
    
    return expected

def correct_for_x_chrom(expected, male_n, female_n):
    """ correct mutations rates for sex-chromosome transmission rates
    
    Args:
        expected: gene-based data frame, containing rates for different mutation
            classes.
         male_n: number of trios with male offspring
         female_n: number of trios with female offspring
    
     Returns:
        a dataframe of mutation rates for genes under different mutation
        classes.
    """
    
    # Calculate the number of transmissions for autosomal, male and female
    # transmissions. The number of transmissions from males is equal to the
    # number of  female probands (since only females receive a chrX from their
    # fathers). Likewise, all offspring receive a chrX from their mothers, so
    # the number of transmissions from females equals the number of probands.
    autosomal = 2 * (male_n + female_n)
    female_transmissions = male_n + female_n
    male_transmissions = female_n
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    alpha = 3.4
    male_factor = 2 / (1 + (1 / alpha))
    female_factor = 2 / (1 + alpha)
    
    # correct the non-PAR chrX genes for fewer transmissions and lower rate
    # (dependent on alpha)
    chrX = expected["chrom"].isin(["X", "chrX"])
    x_factor = ((male_transmissions * male_factor) + (female_transmissions * female_factor)) / autosomal
    x_factor = pandas.Series([x_factor] * len(chrX), index=expected.index)
    x_factor[~chrX] = 1
    
    expected["missense_snv"] *= x_factor
    expected["missense_indel"] *= x_factor
    expected["lof_snv"] *= x_factor
    expected["lof_indel"] *= x_factor
    expected["synonymous_snv"] *= x_factor
    
    return expected

def adjust_indel_rates(expected):
    """ adapt indel rates for lower rate estimate from validated de novos
    
    The indel mutation rates from Samocha et al., Nature Genetics 46:944-950
    assume that the overall indel mutation rate is 1.25-fold greater than the
    overall nonsense mutation rate, ie there are 1.25 times as many frameshifts
    as nonsense mutations. We have our own estimates for the ratio, derived from
    our de novo validation efforts, which we shall apply in place of the Samocha
    et al ratios.
    
    Args:
        rates: data frame of mutation rates.
    
    Returns:
        the rates data frame, with adjusted indel rates.
    """
    
    # the following numbers were derived from the DDD 4K dataset.
    nonsense_n = 411
    frameshift_n = 610
    ddd_ratio = frameshift_n / nonsense_n
    samocha_ratio = 1.25  # Nature Genetics 46:944-950 frameshift to nonsense ratio
    
    # correct back from using the Samocha et al. ratio
    expected["missense_indel"] /= samocha_ratio
    expected["lof_indel"] /= samocha_ratio
    
    # adjust the indel rates for the DDD indel ratio
    expected["missense_indel"] *= ddd_ratio
    expected["lof_indel"] *= ddd_ratio
    
    return expected
