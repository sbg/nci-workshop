# Gene Expression Matrix

A script for merging TCGA Level 3 RNAseq gene expression data (raw counts), along with the cwl.json files for use on the Seven Bridges Cancer Genomics Cloud.

To test locally, use: `python scripts/munger.py -f data/*.txt -c`

For testing the use of an index file: `python scripts/munger.py -r data/test.index -c`