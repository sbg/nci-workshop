# How to use Dockerfiles

# Start with the release_base image from Bioconductor
FROM bioconductor/release_base

MAINTAINER "Gaurav Kaushik" <gaurav@sevenbridges.com>

# Update the container
RUN apt-get update --yes
RUN apt-get install vim --yes

# To add more Bioconductor packages, insert that package in "" in the list below
# More instructions for installing new Bioconductor packages can be found here: https://www.bioconductor.org/install/
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("BiocInstaller", "DESeq2", "limma", "edgeR", "RColorBrewer"))'

# Add the pheatmap package
RUN Rscript -e 'install.packages("pheatmap")'

# Copy the R script and test files
ADD diff.R /usr/local/bin/
ADD brca* ./
RUN chmod a+x /usr/local/bin/diff.R
# When testing the tool, modify diff.R (see comments ~line 18)

# Set CMD to fire up the terminal when you run the container
CMD ["/bin/bash"]
