#' @include dk
#' @author Fabio M D'Orazio
#' @description compile functions

##### FUNCTIONS #####

## extracts tag clusters for all the samples
#iqs <- lapply(sample.names, function(x){
#  tc <- tagClusters(myCAGEset, sample = x, returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
#  tc <- tc[tc$chr %in% paste("chr",1:25, sep = ""),]
#  tc$sampleID = x
#  return(tc)
#})


## annotate cage tag clusters
#' annotate cage tag clusters
#'
#' @param x
#'
#' @return annotated tag clusters
#' @export
#' @importFrom ChIPseeker annotatePeak
#' @examples
annotate.cage.peaks <- function(x){
  cage.range <- toGRanges(x)
  cage.peak <- annotatePeak(cage.range, tssRegion = c(-500, 500), TxDb = txdb,
                            annoDb = 'org.Dr.eg.db', sameStrand = T, verbose = F)
  cage.frame <- data.frame(cage.peak@anno)
  ## get only promoters
  cage.frame.promoters <- subset(cage.frame, cage.frame$annotation == 'Promoter')

  return(cage.frame.promoters)
}

## subset cage object
#' Subset cage object
#'
#' @param x cage.anno.prim
#' @param y up.pgc.high
#' @param z cage.anno.high
#'
#' @return merged data frame of tag clusters
#' @export
#'
#' @examples subset.cage.object(cage.anno.prim, up.pgc.high, cage.anno.high)
subset.cage.object <- function(x,y,z){
  merged.cage.late <- merge(x, y, by.x = 'geneId', by.y = 0)
  merged.cage.early <- merge(z, y, by.x = 'geneId', by.y = 0)

  merged.cage.object <- merge(merged.cage.early, merged.cage.late, by = 'transcriptId')
  merged.cage.object <- subset(merged.cage.object, merged.cage.object$tpm.dominant_ctss.x > 5
                               & merged.cage.object$tpm.dominant_ctss.y > 5)
  ## assign maternal and zygotic dominat peaks
  new.data.frame <- data.frame(seqnames = merged.cage.object$seqnames.x,
                               start = merged.cage.object$start.x, end = merged.cage.object$end.x,
                               strand = merged.cage.object$strand.x,
                               maternal_dominant_ctss = merged.cage.object$dominant_ctss.x,
                               zygotic_dominant_ctss = merged.cage.object$dominant_ctss.y,
                               tpm.maternal_dominant_ctss = merged.cage.object$tpm.dominant_ctss.x,
                               tpm.zygotic_dominant_ctss = merged.cage.object$tpm.dominant_ctss.y
  )
  #sort the clusters by difference between maternal and zygotic dominant peak position
  ##change the sign based on strand
  new.data.frame$difference <- new.data.frame$maternal_dominant_ctss - new.data.frame$zygotic_dominant_ctss
  new.data.frame <- GRanges(new.data.frame)
  new.data.frame$difference[(as.character(strand(new.data.frame)) == '-')] <- -new.data.frame$difference[(as.character(strand(new.data.frame)) == '-')]
  new.data.frame <- data.frame(new.data.frame)
  new.data.frame <- new.data.frame[order(new.data.frame$difference),]

  return(new.data.frame)

}

## PLOT
## creates a new data frame
## centers it either to the maternal dominant or the zygotic dominant
## extract the promoter sequences and plots
#' Title
#'
#' @param merged.cage.object
#' @param method
#'
#' @return
#' @export
#'
#' @examples return.pattern.distr(merged.cage)
return.pattern.distr <- function(merged.cage.object, method = NULL){
  new.data.frame <- data.frame(seqnames = merged.cage.object$seqnames,
                               start = merged.cage.object$start, end = merged.cage.object$end,
                               strand = merged.cage.object$strand,
                               maternal_dominant_ctss = merged.cage.object$dominant_ctss,
                               maternal_IQ_width = merged.cage.object$interquantile_width,
                               zygotic_dominant_ctss = merged.cage.object$dominant_ctss_late,
                               tpm.maternal_dominant_ctss = merged.cage.object$tpm.dominant_ctss,
                               tpm.zygotic_dominant_ctss = merged.cage.object$tpm.dominant_ctss_late,
                               zygotic_IQ_width = merged.cage.object$interquantile_width_late
  )

  ## plot dinucleotide frequncy
  if(method == "maternal") { ## centered to maternal dominant tss
    zebrafishPromotersTSS<-GRanges(seqnames = new.data.frame$seqnames,
                                   ranges=IRanges(start = new.data.frame$maternal_dominant_ctss,
                                                  end = new.data.frame$maternal_dominant_ctss),
                                   strand = new.data.frame$strand,
                                   IQ_width = new.data.frame$maternal_IQ_width,
                                   seqlengths = seqlengths(Drerio))

  }else if(method == "zygotic"){ ## centered to zygotic dominant tss
    zebrafishPromotersTSS<-GRanges(seqnames = new.data.frame$seqnames,
                                   ranges=IRanges(start = new.data.frame$zygotic_dominant_ctss,
                                                  end = new.data.frame$zygotic_dominant_ctss),
                                   strand = new.data.frame$strand,
                                   IQ_width = new.data.frame$zygotic_IQ_width,
                                   seqlengths = seqlengths(Drerio))
  }else{
    stop("'method' parameter must be one of the (\"zygotic\", \"maternal\")")
  }
  ## TATA is between -35 and -23
  ## extend the region upstream of the main CAGE peak with center to -29
  zebrafishPromotersTSSflank <- promoters(zebrafishPromotersTSS, upstream = 50,
                                          downstream = 50)
  ## get 6bp up and 6bp downstream (-35, -23) to include the TATA region
  #zebrafishPromotersTSSflank <- resize(zebrafishPromotersTSSflank, width = 12, fix = 'center')

  return(zebrafishPromotersTSSflank)
}
