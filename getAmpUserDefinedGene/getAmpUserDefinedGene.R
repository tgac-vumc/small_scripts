### load libraries ----
library(QDNAseq)
library(Biobase)

#### PATH AND FILE NAMES ----
input.path <- "user defined input path"
called.file <- file.path(input.path,"30kbp-called.rds" )
gene.file <- file.path(input.path, "target_genes.csv")
output.path <- "user defined output path"

### LOAD DATA ----
data <- readRDS(file = called.file)
geneList <- read.table(file = gene.file, header = TRUE, sep = ";", skip = 0)

### FUNCTIONS ----
.getFeature <- function(calledData, chromosome, position){
  x <- data@featureData@data$chromosome == chromosome & 
             position >= data@featureData@data$start &
             position <= data@featureData@data$end 
  return(x)
}
getGeneCalls <- function(calledData, geneList){
  colnames(geneList) <- tolower(colnames(geneList))
  header <- c("name","chromosome","start")
  if(sum(header %in% colnames(geneList)) != 3){
    stop("GeneList does not contain the required columns.\n  * name\n  * chromosome\n  * start\n\n ")
  }
  nSamples <- length(sampleNames(data))
  results.feature <- matrix(NA, nrow = nrow(geneList), ncol = 3)
  colnames(results.feature) <- c("Chromosome","Start", "End")
  results.calls <- matrix(NA, nrow = nrow(geneList), ncol = nSamples)
  colnames(results.calls) <- sampleNames(data)
  results.probamp <- results.calls
  sel <- rep(FALSE, nrow(geneList))
  for(gene in 1:nrow(geneList)){
    feature <- .getFeature(calledData = data, chromosome = geneList$chromosome[gene], position = geneList$start[gene])
    results.feature[gene,] <- unlist(data@featureData@data[feature,][1:3])
    results.calls[gene,] <- data@assayData$calls[feature,]
    results.probamp[gene,] <- data@assayData$probamp[feature,]
  }
  rownames(results.feature) <- rownames(results.calls) <- rownames(results.probamp) <- geneList$name
  results <- list(sampleNames = sampleNames(data),
                  genes = geneList,
                  features = results.feature,
                  calls = results.calls,
                  probamp = results.probamp)
  return(results)
}
getAmpUserDefinedGene <- function(calledData, geneList, out, col = "green"){
  cat(as.character(Sys.time()), ": Extracting gene data.\n")
  geneData <- getGeneCalls(calledData = data, geneList = geneList)
  for(i in 1:nSamples){
    cat(as.character(Sys.time()), ": Processing sample -", sampleNames(data)[i],"\n")
    sel <- geneData$calls[,i] %in% 2
    sampleData <- cbind(geneData$genes[,1:2],
                        geneData$features,
                        calls = geneData$calls[,i],
                        probamp = geneData$probamp[,i])
    
    output.file <- file.path(out, paste(geneData$sampleNames[i], "_results.csv", sep = ""))
    write.table(x = sampleData, file = output.file, quote = FALSgenE, sep = ";", row.names = FALSE)
    if(TRUE %in% sel == TRUE){
      cat(as.character(Sys.time()), ": Plotting amplifications.\n")
      for(g in 1:nrow(geneData$genes)){
        if(geneData$calls[g,i] == 2){
          output.plot <- paste(geneData$sampleNames[i],"_",
                               geneData$genes[g,1],"_",
                               paste("chr",geneData$features[g,1],"_",geneData$features[g,2],"-",geneData$features[g,3], sep = ""),
                               ".png", sep = "")
          output.plot <- file.path(out, output.plot)
          png(filename = output.plot, width = 1500, height = 750)
          plot(data[chromosomes(data) == geneData$features[g,1],i], gaincol = FALSE, losscol = FALSE, delcol = FALSE, ampcol = col)
          mtext(text = rownames(geneData$features)[g], at = geneData$features[g,2], padj = -1, cex = 0.8)
          dev.off()
        }
      }
    } else {
      cat(as.character(Sys.time()), ": No Amplication calls detected.\n")
    }
  }
}

### RUN THE CODE ----
getAmpUserDefinedGene(calledData = data, geneList = geneList, out = pad, col = "lightblue")

### END ----
