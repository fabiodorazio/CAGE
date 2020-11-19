### Import mapped CAGE bams and converts them to big wigs

library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)

setwd("../newRun/")
files <- list.files(pattern = '\\.sorted.bam$')
## use   coverage.range <- keepSeqlevels(coverage.range, c(1:25), pruning.mode = 'coarse') for zebrafish development samples
save.path <- '~/CAGE_PGC_BED/'
chr <- paste0('chr', c(1:25))

readAlignBAMplus <- function(x){
  
  alignment = readGAlignments(file = x)
  grange <- GRanges(alignment)
  grange.plus <- grange[grange@strand == '+']
  coverage <- coverage(grange.plus)
  coverage.range <- GRanges(coverage)
  coverage.range <- resize(coverage.range, 1, fix = 'start')
  levels(coverage.range@strand) <- c('+', '+', '+')
  
  coverage.range <- keepSeqlevels(coverage.range, chr, pruning.mode = 'coarse')
  #chr <- paste0('chr', c(1:25))
  #seqlevels(coverage.range) <- chr
  
  export.bw(object = coverage.range, con=paste0(save.path, x, '.plus.bw'))
  
}
  

readAlignBAMminus <- function(x){
  
  alignment = readGAlignments(file = x)
  grange <- GRanges(alignment)
  grange.plus <- grange[grange@strand == '-']
  coverage <- coverage(grange.plus)
  coverage.range <- GRanges(coverage)
  coverage.range <- resize(coverage.range, 1, fix = 'start')
  levels(coverage.range@strand) <- c('-', '-', '-')
  
  coverage.range <- keepSeqlevels(coverage.range, chr, pruning.mode = 'coarse')
  seqlevels(coverage.range) <- chr
  
  export.bw(object = coverage.range, con=paste0(save.path, x,'.minus.bw'))
  
}


lapply(files, readAlignBAMplus)
lapply(files, readAlignBAMminus)



