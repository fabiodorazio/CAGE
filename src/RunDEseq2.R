#for pieces of code
## for notes
### titles

### simple function to run DESeq2 given a countdata and input for coldata objects

library(DESeq2)

run.deseq.genes <- function(countdata = countdata, stage, type){
  coldata <- data.frame(stage = stage, condition = type)
  row.names(coldata) <- colnames(totalReadsCount)
  
  dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ stage + condition)
  dds <- DESeq(dds)
  normalizedData <- counts(dds, normalized = TRUE)
  
  ## for testing logFC above or belove a threshold:
  res <- results(dds, lfcThreshold = .5, altHypothesis = 'greaterAbs')
  
  return(res)
}
