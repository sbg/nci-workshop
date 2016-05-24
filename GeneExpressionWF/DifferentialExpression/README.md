##Differential Expression
 
`diff.R` is a script for performing a differential expression analysis using the DESeq2 library in R.

**To run in the command line:** 

`Rscript diff.R brca_gene.csv brca_meta.csv csv_filename pdf_filename <rld | vsd>` 

- `brca_gene.csv` is a matrix of genes vs samples
- `brca_meta.csv` is a matrix of samples vs metadata (e.g. sample_type, gender)
- `*_filename` are the desired prefixes for your output files
- `<rld | vsd>` is your choice of normalization method (use "rld" or "vsd" without quotes)

**Note:** `brca_gene.csv` and `brca_meta.csv` were produced on the Seven Bridges Cancer Genomics Cloud with the "Gene Expression Munger" CWL tool with Level 3 TCGA Gene Quantification Data as the input (see `brca_index_file.txt` the for list of input files). `diff.R` has also been tested with TCGA Exon Quantification data.

###Build your own container
To package this tool into a Docker container for testing and/or deployment, you can use the following command:
`docker build -t <repository/container:tag> .` 

Or you can pull the container using:
`docker pull gauravkaushik/rnaseq:workshop`

###Import into the CGC
You can import **`diff.cwl.json`** in the Seven Bridges CGC to use for your own research or try wrapping this application yourself using the CGC Tool Editor!
