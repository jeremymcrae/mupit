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
import json
import os

from pandas import DataFrame

from mupit.write_json_probands import write_probands_by_gene

class TestWriteJsonProbands(unittest.TestCase):
    
    def setUp(self):
        """ make sure we have a temporary file for writing to
        """
        
        handle, self.path = tempfile.mkstemp()
        self.handle = os.fdopen(handle, 'w')
        
        self.de_novos = DataFrame({
            "hgnc": ["ARID1B", "ARID1B", "AMELX"],
            "person_id": ["id_1", "id_2", "id_1"],
            "consequence": ["missense_variant", "missense_variant",
                "missense_variant"],
            })
    
    def tearDown(self):
        """ ensure we delete temporary files
        """
        
        self.handle.close()
        os.remove(self.path)
    
    def test_write_probands_by_gene(self):
        """ test correct output of json formatted file of probands by gene
        """
        
        write_probands_by_gene(self.de_novos, self.path)
        with open(self.path) as handle:
            self.assertEqual(json.load(handle),
                {'ARID1B': ['id_1', 'id_2'], 'AMELX': ['id_1']})
    
    def test_write_probands_removing_nonfunctional(self):
        """ check that variants that don't alter protein sequence are removed
        """
        
        de_novos = DataFrame({
            "hgnc": ["ARID1B", "ARID1B", "AMELX"],
            "person_id": ["id_1", "id_2", "id_1"],
            "consequence": ["missense_variant", "synonymous_variant",
                "missense_variant"],
            })
        
        write_probands_by_gene(de_novos, self.path)
        with open(self.path) as handle:
            self.assertEqual(json.load(handle),
                {'ARID1B': ['id_1'], 'AMELX': ['id_1']})
    
    def test_write_probands_by_gene_with_handle(self):
        """ check that we can pass a file handle into the json writing function
        """
        
        write_probands_by_gene(self.de_novos, self.handle)
        self.handle.flush()
        
        with open(self.path) as handle:
            self.assertEqual(json.load(handle),
                {'ARID1B': ['id_1', 'id_2'], 'AMELX': ['id_1']})
    
