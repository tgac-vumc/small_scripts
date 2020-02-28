## getAmpUserDefinedGene (Get Amplification of User Define Gene)
### Hendrik F. van Essen (2020)

you need to specify the Path, 'called.Rds' file, and file with Target Genes

install.packages("QDNAseq")
install.packages("Biobase")

## getAmpUserDefinedGene
The *getGeneCalls* function will let the user extract all the Calls from a QDNAseq '.Rds' file with the list of Target Genes as a template.

The input for *getGeneCalls* is the following:
* calledData = QDNAseq called.Rds file
* geneList = data frame with the following columns
** name, hromosome, start

Data will be returned from the function as a list containing the following:
* x$sampleNames : contains the sampleNames as a list
* x$genes : contains the data frame provide by the user
* x$features : contains a data frame of the features from the QDNAseq object where the genes are located
* x$calls : contains a data frame of the calls for each sample for the provided target genes
* x$probamp : contains a data frame of the probablity store for the Amplification Calls for each sample for the provided target genes

The function *getAmpUserDefinedGene* will provide the following for each sample in the data provided.
* .csv file with target genes, calls, probamp
* .png file for each target gene with a positive amplification call (2)

The input for *getAmpUserDefinedGene* is the following:
* calledData = QDNAseq called.Rds file
* geneList = list structure from the function *getGeneCalls*
* out = output path where data needs to go
* col = user define colour of the Amplification calls

### Installation
library(QDNAseq)
library(Biobase)

load the 3 functions into R.
* .getFeature
* getGeneCalls
* getAmpUserDefinedGene

## PIPELINE EXAMPLE
### set paths and file locations
input.path <- "user defined input path"
called.file <- file.path(path,"30kbp-called.rds" )
gene.file <- file.path(pad, "target_genes.csv")
output.path <- "user defined output path"

### load the data 
data <- readRDS(file = called.file)
geneList <- read.table(file = gene.file, header = TRUE, sep = ";", skip = 0)

### run the function
getAmpUserDefinedGene(calledData = data, geneList = geneList, out = output.path, col = "lightblue")

## Citation
Not available
