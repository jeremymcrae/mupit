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

import os
import unittest
import math
import tempfile

import numpy
from pandas import DataFrame, Series

from mupit.gene_enrichment import analyse_enrichment, gene_enrichment, test_enrich
from tests.compare_dataframes import CompareTables

class TestGeneEnrichmentPy(CompareTables):
    """ unit test the gene enrichment testing functions
    """
    
    def setUp(self):
        """ construct some objects for unit tests
        """
        
        self.observed = DataFrame({
            "hgnc": ["ARID1B", "KMT2A"],
            "chrom": ["6", "11"],
            "start_pos": [157150547, 118367083],
            "lof_indel": [1, 0],
            "lof_snv": [1, 0],
            "missense_indel": [1, 0],
            "missense_snv": [1, 1]
            })
        
        self.expected = DataFrame({
            "hgnc": ["ARID1B", "KMT2A"],
            "chrom": ["6", "11"],
            "start_pos": [157150547, 118367083],
            "lof_indel": [0.01, 0.01],
            "lof_snv": [0.005, 0.005],
            "missense_indel": [0.05, 0.05],
            "missense_snv": [0.005, 0.005]
            })
    
    def test_analyse_enrichment(self):
        ''' test that analyse_enrichment works correctly
        '''
        
        de_novos = DataFrame({
            'person_id': ['1'],
            'hgnc': ['KMT2A'],
            'chrom': ['11'], 'start_pos': [1],
            'ref': ['A'], 'alt': ['A'],
            'consequence': ['missense_variant']
            })
        
        rates = DataFrame({
            'hgnc': ['ARID1B', 'KMT2A', 'AMELX'],
            'chrom': ['6', '11', 'X'],
            'syn': [1e-5, 1e-5, 1e-5],
            'mis': [1e-6, 1e-6, 1e-6],
            'non': [1e-7, 1e-7, 1e-7],
            'splice_site': [1e-5, 1e-5, 1e-5],
            'frameshift': [1e-7, 1e-7, 1e-7]
            })
        
        expected = DataFrame({
            'hgnc': ['KMT2A'], 'chrom': ['11'],
            'lof_indel': [0], 'lof_snv': [0],
            'missense_indel': [0], 'missense_snv': [1],
            'p_func': [0.00894529269978444619], 'p_lof': [1.0],
            })
        
        trios = {'male': 200, 'female': 200}
        output = analyse_enrichment(de_novos, trios, rates)
        
        self.compare_tables(output, expected)
        
        # check that we can save a pdf of a Manhattan plot. Save to a temporary
        # file. First check that the file size is zero
        path = tempfile.NamedTemporaryFile()
        self.assertEqual(os.path.getsize(path.name), 0)
        
        # Run the function, then check if the output size is nonzero
        output = analyse_enrichment(de_novos, trios, rates, path)
        self.assertGreater(os.path.getsize(path.name), 0)
    
    def test_lof_enrichment(self):
        """ check the test of loss-of-function enrichment
        """
        
        lof_columns = ["lof_indel", "lof_snv"]
        computed = test_enrich(self.expected, self.observed, lof_columns)
        self.assertEqual(list(computed), [1.113813028913985e-04, 1.00000000])
    
    def test_functional_enrichment(self):
        """ check the test of functional (LoF + missense) enrichment
        """
        
        func_columns = ["lof_indel", "lof_snv", "missense_indel", "missense_snv"]
        computed = test_enrich(self.expected, self.observed, func_columns)
        self.assertEqual(list(computed), [9.4599516119858809e-07, 6.7606180094051796e-02])
    
    def test_gene_enrichment(self):
        """ check that LoF and func P-value columns are added to the dataframe
        """
        
        # add the p-value columns in the expected order.
        expected = self.observed.copy()
        expected["p_lof"] = [1.113813028913985e-04, 1.00000000]
        expected["p_func"] = [9.4599516119858809e-07, 6.7606180094051796e-02]
        
        computed = gene_enrichment(self.expected, self.observed)
        self.compare_tables(computed, expected)
