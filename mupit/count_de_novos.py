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

import numpy
import pandas

from mupit.constants import LOF_CQ, MISSENSE_CQ

def get_de_novo_counts(de_novos):
    """ tallies the mutation types observed for each gene
    
    Args:
        de_novos: data frame listing all the de novo mutations, with columns
            for HGNC symbol, consequence type (VEP style predictions), and a
            column indicating SNV, or indel.
    
    Returns:
        data frame with tally of de novo mutations for each of the mutation
        types.
    """
    
    de_novos = de_novos.copy()
    
    # group the lof and missence consequence strings, and drop all the
    # non-functional de novos
    consequence = de_novos["consequence"].copy()
    consequence[consequence.isin(LOF_CQ)] = "lof"
    consequence[consequence.isin(MISSENSE_CQ)] = "missense"
    de_novos["consequence"] = consequence
    de_novos = de_novos[de_novos["consequence"].isin(["missense", "lof"])]
    
    de_novos["type"] = get_var_type(de_novos)
    
    # count the number of de novos for each type/consequence combination
    if pandas.__version__ < "0.14.0":
        counts = de_novos.pivot_table(values="person_id", rows=["hgnc"],
            cols=["consequence", "type"], aggfunc=len, fill_value=0)
    else:
        counts = de_novos.pivot_table(values="person_id", index=["hgnc"],
            columns=["consequence", "type"], aggfunc=len, fill_value=0)
    
    counts = tidy_count_data(counts, de_novos)
    
    return counts

def get_var_type(de_novos):
    """ figure out whether a de novo mutation is for a snv, or indel
    
    Args:
        de_novos: pandas dataframe of de novo candidate, including "ref" and
            "alt" columns.
    
    Returns:
        pandas Series indicating whether a site is a "snv" or "indel".
    """
    
    if "type" in de_novos.columns:
        var_type = de_novos["type"]
    else:
        var_type = pandas.Series(["snv"] * len(de_novos), index=de_novos.index)
        ref_length = de_novos["ref"].str.len()
        alt_length = de_novos["alt"].str.len()
        
        var_type[(ref_length != 1) | (alt_length != 1)] = "indel"
    
    return var_type

def tidy_count_data(counts, de_novos):
    """ simplify the counts table, since pandas pivot tables are awkward to use.
    
    Args:
        counts: a hierarchically indexed pandas DataFrame of de novo counts
            across different consequence and variant type categories.
        de_novos: pandas dataframe of de novo variants.
    
    Returns:
        count table, with a flattened index, and standardised columns.
    """
    
    counts.columns = [ '_'.join(col).strip() for col in counts.columns.values ]
    
    # Get the genome coordinates of the de novos at the minimum position for
    # each gene. This is so we can make Manhattan plots of significance by
    # genome location.
    counts["hgnc"] = list(counts.index)
    counts["chrom"] = de_novos.groupby("hgnc")["chrom"].agg(min)
    counts["start_pos"] = de_novos.groupby("hgnc")["start_pos"].agg(min)
    
    # make sure all the required count column are present in the table, so that
    # we can later refer to each of the columns
    for col in ["lof_indel", "lof_snv", "missense_indel", "missense_snv"]:
        if col not in counts.columns:
            counts[col] = 0
    
    # remove the hgnc based index, to simplify the table.
    counts = counts.reset_index(drop=True)
    
    # sort the columns to a standard order
    counts = counts[["hgnc", "chrom", "start_pos", "lof_indel", "lof_snv",
        "missense_indel", "missense_snv"]]
    
    return counts
