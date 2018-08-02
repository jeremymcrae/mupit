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

from scipy.stats import poisson

from mupit.count_de_novos import get_de_novo_counts
from mupit.mutation_rates import get_expected_mutations
from mupit.plot_enrichment import plot_enrichment

import pandas

def analyse_enrichment(de_novos, trios, rates=None, plot_path=None):
    """ analyse whether de novo mutations are enriched in genes
    
    Args:
        de_novos: data frame containing all the observed de novos for all the
            genes
        trios: dictionary of male and female proband counts in the population
        plot_path: path to save enrichment plots to, or None
        rates: gene-based mutation rates data frame, or None
    
    Returns:
        data frame containing results from testing for enrichment of de
        in each gene with de novos in it.
    """
    
    observed = get_de_novo_counts(de_novos)
    expected = get_expected_mutations(rates, trios["male"], trios["female"])
    
    # calculate p values for each gene using the mutation rates
    enrichment = gene_enrichment(expected, observed)
    
    # make a manhattan plot of enrichment P values
    if plot_path is not None:
        num_tests = 18500
        plot_enrichment(enrichment, num_tests, plot_path, p_columns=["p_lof", "p_func"])
    
    # remove the position column (which is only used to be able to locate the
    # gene's position on a chromosome on a Manhattan plot).
    del enrichment["start_pos"]
    
    return enrichment

def gene_enrichment(expected, observed):
    """ tests whether genes are enriched with de novo mutations
    
    Args:
        expected: pandas dataframe of expected numbers of mutations per gene,
            given expected mutation rates for each gene.
        observed: pandas data frame with tally of de novo mutations per gene
            for each of the mutation types: lof_snv, lof_indel, missense_snv,
            missense_indel.
    
    Returns:
        pandas data frame with P-values from testing for enrichment.
    """
    
    lof_columns = ["lof_indel", "lof_snv"]
    func_columns = ["lof_indel", "lof_snv", "missense_snv", "missense_indel"]
    
    enriched = observed.copy()
    enriched["p_lof"] = test_enrich(expected, observed, lof_columns)
    enriched["p_func"] = test_enrich(expected, observed, func_columns)
    
    return enriched

def test_enrich(expected, observed, columns):
    """ tests whether genes are enriched with de novo mutations
    
    Args:
        expected: pandas dataframe of expected numbers of mutations per gene,
            given expected mutation rates for each gene.
        observed: pandas data frame with tally of de novo mutations per gene
            for each of the mutation types: lof_snv, lof_indel, missense_snv,
            missense_indel.
        columns: list of columns to use to calculate enrichment within, such as
            the loss-of-function columns ["lof_snv", "lof_indel"].
    
    Returns:
        pandas Series of P-values from testing for enrichment.
    """
    
    # recode the columns in the expected mutations table, so merging the
    # expected and observed datasets doesn't have conflicting column names.
    expected_columns = [ x + "_expected" for x in columns ]
    rename = dict(zip(columns, expected_columns))
    expected = expected.rename(columns=rename)
    
    if 'hgnc' not in observed:
        observed['hgnc'] = observed['symbol']
    
    enriched = observed.merge(expected, how="left", on=["hgnc", "chrom"])
    
    # account for how different pandas versions sum series with only NA
    kwargs = {}
    if pandas.__version__ >= '0.22.0':
        kwargs = {'min_count': 1}
    
    # sum the observed and expected de novo mutations per gene
    observed = enriched[columns].sum(axis=1, **kwargs)
    expected = enriched[expected_columns].sum(axis=1, **kwargs)
    
    # calculate the probability of getting the observed number of de novos,
    # given the expected rate of mutations.
    return poisson.sf(observed - 1, expected)
