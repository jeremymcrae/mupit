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
import math
import tempfile

from pandas import DataFrame

from mupit.combine_analyses import fishersMethod, combine_tests, \
    combine_enrichment_and_clustering
from tests.compare_dataframes import CompareTables

class TestCombineAnalyses(CompareTables):
    
    def test_fishersMethod(self):
        """ check that fishers combined method calculates correctly
        """
        
        p_values = [0.01, 0.001]
        self.assertAlmostEqual(fishersMethod(p_values), 0.0001251292546)

        # correct for single value. This is off by a minuscule fraction due to
        # float division
        p_values = [0.01]
        self.assertAlmostEqual(fishersMethod(p_values), 0.01)

        # correct for list with NA value
        p_values = [0.01, float('nan')]
        self.assertAlmostEqual(fishersMethod(p_values), 0.01)

        # correct for list with NA value
        p_values = [0.01, 0.001, float('nan')]
        self.assertAlmostEqual(fishersMethod(p_values), 0.0001251292546)

        # correct without any p-values
        p_values = [float('nan'), float('nan')]
        self.assertTrue(math.isnan(fishersMethod(p_values)))

        # correct for list with zero
        p_values = [0, 0.01]
        self.assertEqual(fishersMethod(p_values), 0)
    
    def test_combine_enrichment_and_clustering(self):
        """ check that combine_enrichment_and_clustering datasets work correctly
        """
        
        enriched = DataFrame({
            "hgnc": ["GENE1", "GENE2"],
            "p_lof": [0.01, 0.0001],
            "p_func": [0.01, 0.001],
            })
        
        clust = DataFrame({
            "gene_id": ["GENE1", "GENE1", "GENE2", "GENE2"],
            "mutation_category": ["missense", "nonsense", "missense", "nonsense"],
            "probability": [0.1, 0.1, 0.001, 0.001],
            })
        
        expected = DataFrame({
            "hgnc": ["GENE1", "GENE2"],
            "p_lof": [1e-2, 1e-4],
            "p_func": [0.0100, 0.001],
            "p_missense_clust": [0.100, 0.001],
            "p_nonsense_clust": [0.100, 0.001],
            "p_combined": [7.90775527898e-03, 1.48155105579e-05],
            "p_min": [7.90775527898e-03, 1.48155105579e-05],
            })
        
        self.compare_tables(combine_enrichment_and_clustering(enriched, clust), expected)
    
    def test_combine_tests(self):
        """
        """
        
        enriched = DataFrame({
            "hgnc": ["GENE1", "GENE2"],
            "chrom": ["1", "X"],
            "p_lof": [0.01, 0.0001],
            "p_func": [0.01, 0.001],
            })
        
        clustered = DataFrame({
            "gene_id": ["GENE1", "GENE1", "GENE2", "GENE2"],
            "mutation_category": ["missense", "nonsense", "missense", "nonsense"],
            "probability": [0.1, 0.1, 0.001, 0.001],
            })
        
        phenotype = DataFrame({
            "hgnc": ["GENE1", "GENE2", "GENE3"],
            "hpo_similarity_p_value": [0.1, 0.1, 0.001],
            })
        
        # define some temporary files
        meta_clust = tempfile.NamedTemporaryFile()
        meta_enrich = tempfile.NamedTemporaryFile()
        clust = tempfile.NamedTemporaryFile()
        enrich = tempfile.NamedTemporaryFile()
        pheno = tempfile.NamedTemporaryFile()
        
        # write the dataframes to file
        enriched.to_csv(meta_enrich, sep="\t", index=False)
        enriched.to_csv(enrich, sep="\t", index=False)
        clustered.to_csv(meta_clust, sep="\t", index=False)
        clustered.to_csv(clust, sep="\t", index=False)
        phenotype.to_csv(pheno, sep="\t", index=False)
        
        # make sure all the files have finished writing
        meta_clust.flush()
        meta_enrich.flush()
        clust.flush()
        enrich.flush()
        pheno.flush()
        
        computed = combine_tests(meta_clust.name, meta_enrich.name, clust.name,
            enrich.name, pheno.name)
        
        # define the expected output
        expected = DataFrame({
            "chrom": ["1", "X", float('nan')],
            "hgnc": ["GENE1", "GENE2", "GENE3"],
            "meta.p_lof": [0.01, 0.0001, float('nan')],
            "meta.p_func": [0.01, 0.001, float('nan')],
            "meta.p_missense_clust": [0.100, 0.001, float('nan')],
            "meta.p_nonsense_clust": [0.100, 0.001, float('nan')],
            "meta.p_combined": [7.90775527898e-03, 1.48155105579e-05, float('nan')],
            "meta.p_min": [7.90775527898e-03, 1.48155105579e-05, float('nan')],
            "ddd.p_lof": [0.01, 0.0001, float('nan')],
            "ddd.p_func": [0.01, 0.001, float('nan')],
            "ddd.p_missense_clust": [0.100, 0.001, float('nan')],
            "ddd.p_nonsense_clust": [0.100, 0.001, float('nan')],
            "ddd.p_combined": [7.90775527898e-03, 1.48155105579e-05, float('nan')],
            "ddd.p_min": [7.90775527898e-03, 1.48155105579e-05, float('nan')],
            "ddd.hpo_similarity_p_value": [0.1, 0.1, 0.001],
            "p_min": [7.90775527898e-03, 1.48155105579e-05, float('nan')],
            })
        
        self.compare_tables(computed, expected)
        
