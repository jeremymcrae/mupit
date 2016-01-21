"""
Copyright (c) 2015 Genome Research Ltd.
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

from pandas import DataFrame, Series, pivot_table, MultiIndex

from mupit.count_de_novos import get_de_novo_counts, get_var_type, \
    tidy_count_data

class TestCountDeNovosPy(unittest.TestCase):
    """ unit test the de novo mutation counting functions
    """
    
    def setUp(self):
        """ construct some objects for unit tests
        """
        
        self.de_novos = DataFrame({
            "person_stable_id": ["1"],
            "hgnc": ["KMT2A"],
            "chrom": ["1"],
            "start_pos": [1],
            "ref": ["A"],
            "alt": ["A"],
            "consequence": ["missense_variant"]
            })
        
        self.expected = DataFrame({
            "hgnc": ["KMT2A"],
            "chrom": ["1"],
            "start_pos": [1],
            "lof_indel": [0],
            "lof_snv": [0],
            "missense_indel": [0],
            "missense_snv": [0]
            })
    
    def compare_tables(self, first, second):
        """ compare two pandas Dataframes to see if they contain the same Dataframes
        
        Note that the dataframes don't have to have the same column order.
        
        Args:
            first: pandas DataFrame
            second: pandas DataFrame
        """
        
        # make sure the two tables have the same columns
        self.assertEqual(set(first.columns), set(second.columns))
        
        # make sure the values in the columns are identical between tables
        for column in first:
            self.assertEqual(list(first[column]), list(second[column]))
    
    def test_get_de_novo_counts_missense(self):
        """ check the de novo counts for missense-equivalent consequences
        """
        
        # check all the possible missense-equivalent SNV VEP consequences
        snv = ["missense_variant", "stop_lost", "coding_sequence_variant",
            "protein_altering_variant"]
        
        self.expected["missense_snv"] = 1
        for cq in snv:
            self.de_novos["consequence"] = cq
            computed = get_de_novo_counts(self.de_novos)
            self.compare_tables(computed, self.expected)
        
        # check all the possible missense-equivalent inframe insertion and
        # deletion VEP consequences
        indel = ["inframe_deletion", "inframe_insertion"]
        
        self.expected["missense_snv"] = 0
        self.expected["missense_indel"] = 1
        self.de_novos["ref"] = "GG"
        for cq in indel:
            self.de_novos["consequence"] = cq
            computed = get_de_novo_counts(self.de_novos)
            self.compare_tables(computed, self.expected)
    
    def test_get_de_novo_counts_loss_of_function(self):
        """ check the de novo counts for loss-of-function consequences
        """
        
        # check all the loss-of-function-equivalent SNV VEP consequences
        snv = ["stop_gained", "splice_acceptor_variant", "splice_donor_variant",
             "initiator_codon_variant", "start_lost",
            "conserved_exon_terminus_variant"]
        
        self.expected["lof_snv"] = 1
        for cq in snv:
            self.de_novos["consequence"] = cq
            computed = get_de_novo_counts(self.de_novos)
            self.compare_tables(computed, self.expected)
        
        # check all the loss-of-function-equivalent inframe insertion and
        # deletion VEP consequences
        indel = ["frameshift_variant"]
        
        self.expected["lof_snv"] = 0
        self.expected["lof_indel"] = 1
        self.de_novos["ref"] = "GG"
        for cq in indel:
            self.de_novos["consequence"] = cq
            computed = get_de_novo_counts(self.de_novos)
            self.compare_tables(computed, self.expected)
    
    def test_get_de_novo_counts_nonfunctional(self):
        """ check the de novo counts for nonfunctional consequence types
        """
        
        snv = ["transcript_ablation", "transcript_amplification", \
            "incomplete_terminal_codon_variant", "stop_retained_variant", \
            "synonymous_variant", "mature_miRNA_variant", \
            "5_prime_UTR_variant", "3_prime_UTR_variant", \
            "non_coding_transcript_exon_variant", "intron_variant", \
            "NMD_transcript_variant", "non_coding_transcript_variant", \
            "upstream_gene_variant", "downstream_gene_variant", \
            "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", \
            "regulatory_region_ablation", "regulatory_region_amplification", \
            "feature_elongation", "regulatory_region_variant", \
            "feature_truncation", "intergenic_variant"]
        
        # a dataframe with only nonfunctional consequences should raise an
        # error, due to difficulties in counting rows
        for cq in snv:
            self.de_novos["consequence"] = cq
            with self.assertRaises(KeyError):
                computed = get_de_novo_counts(self.de_novos)
        
        # add a functional variant to the dataframe, so that the counting
        # doesn't raise an error
        new_row = self.de_novos.copy()
        new_row["consequence"] = "missense_variant"
        self.de_novos = self.de_novos.append(new_row, ignore_index=True)
        self.expected["missense_snv"] = 1
        
        for cq in snv:
            self.de_novos["consequence"][0] = cq
            computed = get_de_novo_counts(self.de_novos)
            self.compare_tables(computed, self.expected)
    
    def test_get_de_novo_counts_multiple(self):
        """ check the de novo counts for multiple variants
        """
        
        # add another variant to the dataframe
        new_row = self.de_novos.copy()
        self.de_novos = self.de_novos.append(new_row, ignore_index=True)
        
        # counting two of the same variant give a count of 2
        self.expected["missense_snv"] = 2
        self.de_novos["consequence"] = "missense_variant"
        computed = get_de_novo_counts(self.de_novos)
        self.compare_tables(computed, self.expected)
        
        # counting two different variant types gives two counts of 1
        self.expected["missense_snv"] = 1
        self.expected["lof_snv"] = 1
        self.de_novos["consequence"][0] = "stop_gained"
        computed = get_de_novo_counts(self.de_novos)
        self.compare_tables(computed, self.expected)
        
        # now check when we have two variants in different genes
        new_row = self.expected.copy()
        self.expected = self.expected.append(new_row, ignore_index=True)
        self.de_novos["hgnc"][0] = "ARID1B"
        self.expected["hgnc"][0] = "ARID1B"
        self.expected["missense_snv"][0] = 0
        self.expected["lof_snv"][1] = 0
        computed = get_de_novo_counts(self.de_novos)
        self.compare_tables(computed, self.expected)
    
    def test_get_var_type(self):
        """ check that get_var_type works correctly
        """
        
        # check that if a "type" column is present, it just returns that
        self.de_novos["type"] = "TEST"
        self.assertEqual(list(get_var_type(self.de_novos)), ["TEST"])
        
        del self.de_novos["type"]
        
        # add another variant to the dataframe, then check that sites with
        # various ref lengths return the correct type annotations
        new_row = self.de_novos.copy()
        self.de_novos = self.de_novos.append(new_row, ignore_index=True)
        self.de_novos["ref"] = ["A", "AG"]
        self.assertEqual(list(get_var_type(self.de_novos)), ["snv", "indel"])
        
        # set the ref length so that they could be SNVs, then check that sites
        # with various ALT lengths return the correct type annotations
        self.de_novos["ref"] = ["A", "A"]
        self.de_novos["alt"] = ["AG", "A"]
        self.assertEqual(list(get_var_type(self.de_novos)), ["indel", "snv"])
    
    def test_tidy_count_data(self):
        """ make sure that standardising the count table works correctly
        """
        
        arrays = [[1, 1, 2, 2], ['red', 'blue', 'red', 'blue']]
        index = MultiIndex.from_arrays(arrays, names=['first', 'second'])
        
        
        # counts = DataFrame({})
        
        
