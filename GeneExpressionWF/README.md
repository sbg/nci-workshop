##GeneExpressionWF
This workflow 1) merges any number of Level 3 TCGA RNASeq Gene Expression Quantification files into a single matrix and also produces a "metadata matrix" with information about individual cases (*GeneExpressionMunger*) using Python and 2) runs a custom differential expression analysis that produces a matrix of (adj) p-values per gene based on differences among sample-types (e.g. Solid Tissue Normal and Primary Tumor) as well as some plots in PDF (*DifferentialExpression)*. 

###Import into the CGC
You can import **`gene2diff_workflow.cwl.json`** in the Seven Bridges CGC to use for your own research.