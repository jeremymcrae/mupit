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

import unittest

import numpy
from pandas import DataFrame, Series, pivot_table, MultiIndex

from mupit.mutation_rates import get_expected_mutations, correct_for_x_chrom, \
    adjust_indel_rates
from tests.compare_dataframes import CompareTables

class TestMutationRates(CompareTables):
    
    def test_get_expected_mutations(self):
        """ check  that we calculate the expected mutation rates correctly
        """
        
        male = 1000
        female = 1000
        
        rates = DataFrame({
            "hgnc": ["ARID1B", "KMT2A", "AMELX"],
            "chrom": ["6", "11", "X"],
            "syn": [1e-5, 1e-5, 1e-5],
            "mis": [1e-6, 1e-6, 1e-6],
            "non": [1e-7, 1e-7, 1e-7],
            "splice_site": [1e-5, 1e-5, 1e-5],
            "frameshift": [1e-7, 1e-7, 1e-7]
            })
        
        expected = DataFrame({
            "hgnc": ["ARID1B", "KMT2A", "AMELX"],
            "chrom": ["6", "11", "X"],
            "missense_snv": [0.004000000000, 0.004000000000, 0.0024545454545454545],
            "lof_snv": [0.040400000000000005, 0.040400000000000005, 0.024790909090909096],
            "missense_indel": [0.00005004887585532747, 0.00005004887585532747, 0.00003071181018395095],
            "lof_indel": [0.00045043988269794719, 0.00045043988269794719, 0.00027640629165555852],
            "synonymous_snv": [0.040000, 0.040000, 0.024545454545454547],
            })
    
        
        computed = get_expected_mutations(rates, male, female)
        self.compare_tables(computed, expected)
    
    def test_correct_for_x_chrom(self):
        """ check that correcting for differences in chrX mutation rates is correct
        """
        
        rates = DataFrame({
            "hgnc": ["ARID1B", "KMT2A", "AMELX"],
            "chrom": ["6", "11", "X"],
            "missense_snv": [0.001, 0.001, 0.001],
            "lof_snv": [0.001, 0.001, 0.001],
            "missense_indel": [0.001, 0.001, 0.001],
            "lof_indel": [0.001, 0.001, 0.001],
            "synonymous_snv": [0.001, 0.001, 0.001],
            })
        
        expected = DataFrame({
            "hgnc": ["ARID1B", "KMT2A", "AMELX"],
            "chrom": ["6", "11", "X"],
            "missense_snv": [0.001, 0.001, 0.00061363636363636362],
            "lof_snv": [0.001, 0.001, 0.00061363636363636362],
            "missense_indel": [0.001, 0.001, 0.00061363636363636362],
            "lof_indel": [0.001, 0.001, 0.00061363636363636362],
            "synonymous_snv": [0.001, 0.001, 0.00061363636363636362],
            })
        
        # check with even balance of males and females
        self.compare_tables(correct_for_x_chrom(rates, 1000, 1000), expected)
        
        # check that rates are unchanged if no males are present
        self.compare_tables(correct_for_x_chrom(rates, 0, 1000), rates)
    
    def test_adjust_indel_rates(self):
        """ check that adjusting indel rates gives correct results
        """
        
        rates = DataFrame({
            "hgnc": ["ARID1B", "KMT2A", "AMELX"],
            "chrom": ["6", "11", "X"],
            "missense_snv": [0.001, 0.001, 0.001],
            "lof_snv": [0.001, 0.001, 0.001],
            "missense_indel": [0.001, 0.001, 0.001],
            "lof_indel": [0.001, 0.001, 0.001],
            "synonymous_snv": [0.001, 0.001, 0.001],
            })
        
        expected = DataFrame({
            "hgnc": ["ARID1B", "KMT2A", "AMELX"],
            "chrom": ["6", "11", "X"],
            "missense_snv": [0.001, 0.001, 0.001],
            "lof_snv": [0.001, 0.001, 0.001],
            "missense_indel": [0.0011260997067448681, 0.0011260997067448681, 0.0011260997067448681],
            "lof_indel": [0.0011260997067448681, 0.0011260997067448681, 0.0011260997067448681],
            "synonymous_snv": [0.001, 0.001, 0.001],
            })
        
        self.compare_tables(adjust_indel_rates(rates), expected)
