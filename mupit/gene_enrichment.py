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
        columns: list of ciolumns to use to calculate enrichment within, such as
            the loss-of-fuinction columns ["lof_snv", "lof_indel"].
    
    Returns:
        pandas Series of P-values from testing for enrichment.
    """
    
    # recode the columns in the expected mutations table, so merging the
    # expected and observed datasets doesn't have conflicting column names.
    expected_columns = [ x + "_expected" for x in columns ]
    rename = dict(zip(columns, expected_columns))
    expected = expected.rename(columns=rename)
    
    enriched = observed.merge(expected, how="left", on=["hgnc", "chrom"])
    
    # sum the observed and expected de novo mutations per gene
    observed = enriched[columns].sum(axis=1)
    expected = enriched[expected_columns].sum(axis=1)
    
    # calculate the probability of getting the observed number of de novos,
    # given the expected rate of mutations.
    return poisson.sf(observed - 1, expected)
