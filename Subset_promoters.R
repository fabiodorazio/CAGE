## subset sharp or broad promoters
.extract.sharp.promoters <- function(x, method = 'sharp', width = 100, fix = Fix){
  sub.promoters <- x
    if(method == "sharp") { ## centered to maternal dominant tss
        sharp.promoters <- subset(sub.promoters, sub.promoters$interquantile_width <= 9)
  }else if(method == "broad"){ ## centered to zygotic dominant tss
                sharp.promoters <- subset(sub.promoters, sub.promoters$interquantile_width > 9)
  }else{
    stop("'method' parameter must be one of the (\"sharp\", \"broad\")")
  }

      sharp.promoters.range <- GRanges(seqnames = sharp.promoters$chr, IRanges(start = sharp.promoters$dominant_ctss, end = sharp.promoters$dominant_ctss), strand = sharp.promoters$strand)
  sharp.promoters.range <- resize(sharp.promoters.range, width, fix = fix)
  return(sharp.promoters.range)
}
