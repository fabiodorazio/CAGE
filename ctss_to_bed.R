## convert .ctss file in .bed

ctss.to.bed <- function(x){
  
  file <- read.csv(x, sep = '\t', header = F)
  
  file.sort <- file[,c(1,3,6,5)]
  colnames(file.sort) <- c('chr', 'pos', 'strand', 'count')
  file.sort$pos <- file.sort[,2] + 1
  write.table(file.sort, paste0(x, 'converted.bed'), quote=F, sep="\t", row.names=F, col.names=F)
}
file.bed = list.files(pattern = '\\.ctss.bed')
All.bed.file = lapply(file.bed, ctss.to.bed)
