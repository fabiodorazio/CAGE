## subset cage sample

subset.cage.sample <- function(x,y){
  x <- resize(x,fix='center', width=width(x)+100) ## 50 up and 50 downstream
  y <- resize(y, fix='center', width=width(y)+100)
  ov.late <- findOverlaps(x, y) 
  x[queryHits(ov.late),]
  ov.late <- y[subjectHits(ov.late),]
  ov.late <- data.frame(ov.late)
  colnames(ov.late) <- paste0(colnames(ov.late), '_late')
  
  ov.early <- findOverlaps(y, x) 
  y[queryHits(ov.early),]
  ov.early <- x[subjectHits(ov.early),]
  ov.early <- data.frame(ov.early)
  ## do not rename the columns to allow granges to find the coordinates
  
  merged.late.early.frame <- cbind(ov.early, ov.late)
  merged.late.early <- GRanges(merged.late.early.frame)
  
  ## annotate peaks
  ## requires cage.pipe function
  source('../cage_pipe.R')
  cage.anno <- cage.pipe.tbp2(merged.late.early, 'org.Dr.eg.db')
  #cage.anno <- subset(cage.anno, cage.anno$tpm.dominant_ctss > 3 & cage.anno$tpm.dominant_ctss_late > 3)
  cage.anno$difference <- cage.anno$dominant_ctss_late - cage.anno$dominant_ctss
  ## change the sign based on strand
  cage.anno <- GRanges(cage.anno)
  cage.anno$difference[(as.character(strand(cage.anno)) == '-')] <- -cage.anno$difference[(as.character(strand(cage.anno)) == '-')]
  cage.anno <- data.frame(cage.anno)
  
  return(cage.anno)
}

