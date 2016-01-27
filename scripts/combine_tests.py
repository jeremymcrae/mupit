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

from mupit.combine_analyses import combine_tests
from mupit.open_ddd_data import open_known_genes

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("--known-genes", help="Path to known_gene file.")
    parser.add_argument("--ddd-enrichment",
        help="Path to results from enrichment testing of DDD de novos.")
    parser.add_argument("--meta-enrichment",
        help="Path to results from enrichment testing of meta-analysis de novos.")
    parser.add_argument("--ddd-clustering",
        help="Path to results from proximity clustering testing of DDD de novos.")
    parser.add_argument("--meta-clustering",
        help="Path to results from proximity clustering testing of meta-analysis de novos.")
    parser.add_argument("--ddd-phenotype",
        help="Path to file of DDD HPO similarity testing results.")
    parser.add_argument("--output",
        help="Path to write output file for combined testing results.")
    
    args = parser.parse_args()
    
    return args

def include_known_gene_status(merged, known_gene_path):
    """ annotate genes with their known developmental disorder status
    """
    
    # load the table of known developmental disorder genes
    known_gene = open_known_genes(known_gene_path)
    
    # annotate each column with DDG2P status
    merged["known"] = merged["hgnc"].isin(known_gene["gene"])
    merged["known_dominant"] = merged["hgnc"].isin(known_gene["gene"][known_gene["dominant"]])
    
    return merged

def main():
    args = get_options()
    
    merged = combine_tests(args.meta_clustering, args.meta_enrichment,
        args.ddd_clustering, args.ddd_enrichment, args.ddd_phenotype)
    merged = include_known_gene_status(merged, args.known_genes)
    merged.to_csv(args.output, sep="\t", index=False, na_rep="NA")

if __name__ == '__main__':
    main()
