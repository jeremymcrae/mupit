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
import tempfile

from pandas import DataFrame

from mupit.open_ddd_data import standardise_ddd_de_novos, open_known_genes
from tests.compare_dataframes import CompareTables

class TestOpenDddData(CompareTables):
    
    def test_standardise_ddd_de_novos(self):
        """ test correct output of json formatted file of probands by gene
        """
        
        temp = tempfile.NamedTemporaryFile()
        variants = DataFrame({'person_stable_id': ['id1', 'id2'],
            'sex': ['F', 'M'],
            'chrom': ['1', 'X'],
            'pos': [1000, 2000],
            'ref': ['A', 'G'],
            'alt': ['T', 'GGG'],
            'symbol': ['TEMP1', 'TEMP2'],
            'consequence': ['missense_variant', 'frameshift_variant'],
        })
        
        expected = DataFrame({'person_id': ['id1', 'id2'],
            'sex': ['female', 'male'],
            'chrom': ['1', 'X'],
            'start_pos': [1000, 2000],
            'end_pos': [1000, 2000],
            'ref_allele': ['A', 'G'],
            'alt_allele': ['T', 'GGG'],
            'hgnc': ['TEMP1', 'TEMP2'],
            'consequence': ['missense_variant', 'frameshift_variant'],
            'type': ['snv', 'indel'],
            'study_code': ['ddd_unpublished', 'ddd_unpublished'],
            'publication_doi': [None, None],
            'study_phenotype': ['developmental_disorders', 'developmental_disorders'],
        })
        
        variants.to_csv(temp.name, sep='\t', index=False)
        observed = standardise_ddd_de_novos(temp.name)
        self.compare_tables(observed, expected)
    
    def test_open_known_genes(self):
        """ check that opening known DD genes works correctly
        """
        
        temp = tempfile.NamedTemporaryFile()
        genes = DataFrame({
            'gene': ['TEMP1', 'TEMP2'],
            'ddg2p_status': ['Confirmed DD Gene', 'Possible DD Gene'],
            'mode': ['Monoallelic', 'X-linked dominant'],
            'mech': ['loss-of-function', 'Activating'],
            'chr': ['chr1', '1'],
            'start': [1000, 2000],
            'end': [1500, 2500],
        })
        
        expected = DataFrame({
            'gene': ['TEMP1'],
            'type': ['Confirmed DD Gene'],
            'mode': ['Monoallelic'],
            'mech': ['Loss-of-function'],
            'chr': ['1'],
            'start': [1000],
            'end': [1500],
            'dominant': [True],
            'hemizygous': [False],
        })
        
        genes.to_csv(temp.name, sep='\t', index=False)
        observed = open_known_genes(temp.name)
        self.compare_tables(observed, expected)
        
        # and check that it still works cleanly with the pipe separator
        genes.to_csv(temp.name, sep='|', index=False)
        observed = open_known_genes(temp.name)
        self.compare_tables(observed, expected)
    
    def test_open_known_genes_alternate(self):
        """ check that opening known DD genes works correctly with alternate columns
        """
        
        temp = tempfile.NamedTemporaryFile()
        genes = DataFrame({
            'gencode_gene_name': ['TEMP1', 'TEMP2'],
            'type': ['Confirmed DD Gene', 'Possible DD Gene'],
            'mode': ['Monoallelic', 'X-linked dominant'],
            'mech': ['loss-of-function', 'Activating'],
            'chr': ['chr1', '1'],
            'start': [1000, 2000],
            'end': [1500, 2500],
        })
        
        expected = DataFrame({
            'gene': ['TEMP1'],
            'type': ['Confirmed DD Gene'],
            'mode': ['Monoallelic'],
            'mech': ['Loss-of-function'],
            'chr': ['1'],
            'start': [1000],
            'end': [1500],
            'dominant': [True],
            'hemizygous': [False],
        })
        
        genes.to_csv(temp.name, sep='\t', index=False)
        observed = open_known_genes(temp.name)
        self.compare_tables(observed, expected)
    
    
