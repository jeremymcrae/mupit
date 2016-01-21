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

import json

def write_probands_by_gene(de_novos, path):
    """ Write a json-formatted list of probands per gene
    
    We want a correctly-formatted JSON list of probands per gene, for analysis
    of phenotype similarity between probands, see
    https://github.com/jeremymcrae/hpo_similarity.
    
    Args:
        de_novos: dataframe of de-novo mutations per proband
        path: path to write the json data to.
    """
    
    probands_by_gene = {}
    for (gene, group) in de_novos.groupby("hgnc", sort=True):
        probands_by_gene[gene] = sorted(group["person_id"])
    
    with open(path) as handle:
        json.dumps(probands_by_gene, handle, indent=True)
