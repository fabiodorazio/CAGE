### merge two consensus clusters based on granges overlap

myCAGEsetConsensus <- readRDS('CAGEset_PGC_soma_Early_Late.rds')
aggregateTagClusters(myCAGEsetConsensus)

sample.names <- unname(sampleLabels(myCAGEsetConsensus))
iqs <- lapply(sample.names, function(x){
  tc <- tagClusters(myCAGEsetConsensus, sample = x, returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
  tc <- tc[tc$chr %in% paste("chr",1:25, sep = ""),]
  tc$sampleID = x
  return(tc)
})

iqs_PGC <- iqs[c(1,3)] ## subset only early and late PGCs
iqs_earlyPGC <- GRanges(as.data.frame(iqs_PGC[1]))
iqs_latePGC <- GRanges(as.data.frame(iqs_PGC[2]))

iqs_gg <- bind_rows(iqs)

tc_gr <- makeGRangesFromDataFrame(iqs_gg, keep.extra.columns = TRUE)
## subset promoter peaks (merge it with promoters_df)
promoters.range <- GRanges(promoters_df)
iqs_gg.range <- GRanges(iqs_gg)
iqs_gg_promoters <- subsetByOverlaps(iqs_gg.range, promoters.range)

subset_granges_for_shift <- function(x,y){
  x <- resize(x,fix='center', width=width(x)+10)
  y <- resize(y, fix='center', width=width(y)+10)
  ov.late <- findOverlaps(x, y) ## finds overlaps and link each query in grange1 to each in grange2 
  x[queryHits(ov.late),] 
  ov.late <- y[subjectHits(ov.late),]
  ov.late <- data.frame(ov.late)
  colnames(ov.late) <- paste0(colnames(ov.late), '_late')
  
  ov.early <- findOverlaps(y, x) 
  y[queryHits(ov.early),]
  ov.early <- x[subjectHits(ov.early),]
  ov.early <- data.frame(ov.early)
  
  merged.late.early.frame <- cbind(ov.early, ov.late)
  merged.late.early.frame <- subset(merged.late.early.frame,
                                    merged.late.early.frame$tpm.dominant_ctss_late > 1 | 
                                      merged.late.early.frame$tpm.dominant_ctss > 1)
  merged.late.early <- GRanges(merged.late.early.frame)
  
}

cage.PGC.merged <- subset_granges_for_shift(iqs_earlyPGC, iqs_latePGC)
