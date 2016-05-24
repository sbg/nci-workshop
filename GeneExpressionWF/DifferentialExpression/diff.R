####
# diff.R is a script that runs a differential expression analysis on RNASeq data (counts)
# To run in the command line:
# Rscript diff.R brca_gene.csv brca_meta.csv csv_filename pdf_filename <rld | vsd>
#   where:  brca_gene.csv is a matrix of genes vs samples
#           brca_meta.csv is a matrix of samples vs metadata (e.g. sample_type, gender)
#           *_filename is the desired prefix for your output files
#           <rld | vsd> is your choice of normalization method (use "rld" or "vsd" without quotes)
###

# Import libraries
library(limma)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
options(warn=-1)
print("Libraries Imported")

# Read CSV files (harcoded for now)
# Collect arguments from the command line
args <- commandArgs(TRUE)
gene <- read.csv(args[1], header=TRUE, row.names=1) # get the counts per gene (row) per sample/case (col)
meta <- read.csv(args[2], header=TRUE) # get the metadata matrix

## For testing only:
#gene = read.csv("brca_gene.csv", header = TRUE, row.names=1)
#meta = read.csv("brca_metadata.csv", header = TRUE)
#gene = gene[,1:20]
#meta = meta[1:20,]

print("CSV Files Read")

# Get the Differential Expression results
print("Initializing Differential Expression Analysis -- go grab a coffee :)")
meta.df = data.frame(meta)
dds = DESeqDataSetFromMatrix(countData = gene, colData = meta.df, ~ sample_type + X.case_id)
design(dds) <- ~ sample_type + X.case_id
dds = DESeq(dds)
#resultsNames(dds)
res <- results(dds, contrast = c("sample_type", "Solid Tissue Normal", "Primary Tumor"))
resOrdered <- res[order(res$padj),] # order the results by adjusted p-value
print("Differential Expression analysis DONE!")

# Output your data report
'%&%' <- function(x, y)paste0(x,y) # create string concat func
output_title <- args[3] %&% "_rnaseq.csv"
write.csv(resOrdered, output_title) # output df as csv

#### PLOTS ####

# Create PDF with custom title
plot_title <- args[4] %&% "_plots_rnaseq.pdf"
pdf(plot_title)

# Plot -- MA
plotMA(res, main="DESeq2", ylim = c(-4,4))

# Set up normalization method
normalization_method <- args[5]
if (normalization_method == "rld") {
    nmeth <- rlog(dds,blind=FALSE)
} else if (normalization_method == "vsd") {
    nmeth <- varianceStabilizingTransformation(dds, blind=FALSE)
}
# PCA plot -- sample-type
plotPCA(nmeth,intgroup=c("sample_type"))

# Finish up
dev.off()
print("PDF plot saved")

######################