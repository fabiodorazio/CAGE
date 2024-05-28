#### SUBSET TSS BASED ON DIFF GENE EXPRESSION ####
## execution:
## - add desired DEG files to the files list and run the code until line 88
## - substitute the variable list name to be analyzed on line 102
## - replace the indexing operator on line matching the desired cage set 129

library(AnnotationDbi)
library(purrr)
library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicFeatures)
library(dplyr)
library(BSgenome.Drerio.UCSC.danRer7)

danRer7CAGEset <- readRDS("CAGEsetPGCsomaEarlyLateMerged.replicates.rsd")


files <- c('DownregulatedInSomaOnly_Prim5vs10somites.txt', 
           'DownregulatedInPGConly_Prim5vs10somites.txt',
           'UpregulatedGenesInPGCsvsSomaAtHigh.txt', 'UpregulatedGenesInPGCvsSomaPrim5.txt',
           'UpregulatedInPGConly_Prim5vs10somites.txt',
           'UpregulatedInSomavsPGCatPrim5.txt',
           'UpregulatedInBothPGCandSoma_Prim5vs10somites.txt',
           'UPREGinPGCfromHighToDome.txt',
           'UpregulatedInSomaAtPrim5vsPGCandSomitesStage.txt',
           'UpregulatedInSomaAtPrim5vsPGCandDomeStageDanRer7-10.txt',
           'upregulatedInSomaAtDomeVSpgcs.txt',
           'UpregulatedInPGCsVSsomaAtprim5AndFromDome.txt',
           'UpregulatedinPGCsFromDomeToPrim5.txt')

.prepare_data <- function(file){
  x <- read.csv(file, sep = '\t')
  x_subset <- subset(x, abs(x$log2FoldChange) > 2)
  if('X' %in% colnames(x)){
    rownames(x_subset) <- as.character(x_subset$X)
    x_subset$X <- NULL
  }
  else{
    x_subset <- x_subset
  }
  return(x_subset)
}
all_files_rna <- lapply(files, .prepare_data)
#name elements
for(element in 1:length(all_files_rna)){
  names(all_files_rna)[element] <- strsplit(files[element], '\\.')[[1]][1]
}

## load Tx data
txdb <-loadDb('~/Desktop/Postdoc/Data_Analysis/annotation/txdb_DanRer7.sqlite')
txdb7 <- toGRanges(txdb)

## merge deseq table with gene coordinates
## matches DGE class and logFoldChange with gene coordinates
merge.with.coordinates <- function(x,y){
  deseq.coordinates <- merge(as.data.frame(txdb7), x, by = 0) %>%
    mutate(cc_id = paste(seqnames,start,end,width,strand,sep = "_"))
  deseq.res.gr <- GRanges(seqnames = deseq.coordinates$seqnames,
                          ranges = IRanges(start = deseq.coordinates$start,
                                           end = deseq.coordinates$end),
                          strand = deseq.coordinates$strand)
  values(deseq.res.gr) <- deseq.coordinates[,1:ncol(deseq.coordinates)]
  
  # extract tss
  a <- tagClusters(danRer7CAGEset, sample = y, returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
  a <- a[a$chr %in% paste("chr",1:25, sep = ""),]
  a$sampleID <- y
  a <- mutate(a, promoterID = paste(chr,start,end,strand, sep = "_"))
  a <- annotatePeak(GRanges(a), TxDb=txdb, tssRegion=c(-1000, 100), verbose=FALSE, 
                    sameStrand = TRUE)
  gr <- a@anno
  gr <- subset(gr, gr$annotation == 'Promoter')
  
  merge.rnaseq.tss <- function(x, y){
    # determine all overlaps and have granges map to TC coordinates
    hits <- findOverlaps(x,y)
    x_sub <- x[queryHits(hits)]
    y_sub <- y[subjectHits(hits)]
    tc_deseq <- x_sub
    values(tc_deseq) <- cbind(values(tc_deseq), values(y_sub))
    # keep only highest:
    df <- tibble::as_tibble(data.frame(promoterID = tc_deseq$promoterID,
                            tpm = tc_deseq$tpm,
                            cc_id = tc_deseq$cc_id)) %>%
      group_by(cc_id) %>%
      dplyr::slice(which.max(tpm))
    # select only these promoterIDs
    sel <- which(tc_deseq$promoterID %in% df$promoterID & tc_deseq$cc_id %in% df$cc_id)	
    return(tc_deseq[sel])
  }
  tcs_res <- merge.rnaseq.tss(x = gr, y = deseq.res.gr)

  return(tcs_res)
}

if (!interactive()){
    PGC_high_CAGE <- lapply(all_files_rna, merge.with.coordinates, y = unname(danRer7CAGEset@sampleLabels)[1])
    Soma_high_CAGE <- lapply(all_files_rna, merge.with.coordinates, y = unname(danRer7CAGEset@sampleLabels)[2])
    PGC_prim5_CAGE <- lapply(all_files_rna, merge.with.coordinates, y = unname(danRer7CAGEset@sampleLabels)[3])
    Soma_prim5_CAGE <- lapply(all_files_rna, merge.with.coordinates, y = unname(danRer7CAGEset@sampleLabels)[4])
}


## order by iq_width 
order_by_iqwidt <- function(x){
  ordered.tc <- x[order(x$interquantile_width),]
  return(ordered.tc)
}
## 
my.deg.list <- lapply(PGC_prim5_CAGE, order_by_iqwidt)


# plot
.retrieveSequences <- function(x, range = c(400,400), SeqLengths, bs_genome, remove.na = FALSE){
  x <- data.frame(x)
  x <- x[sample(1:nrow(x), nrow(x)),] # rewrite x by shuffle sequences 
  # into granges from dominant
  tc.gr <- GRanges(seqnames = x$seqnames, ranges = IRanges(start = x$dominant_ctss, end = x$dominant_ctss), 
                   strand = x$strand, interquantile_width = x$interquantile_width, seqlengths = SeqLengths)
  # extend grange to the flanking regions
  tc.flank <- trim(promoters(tc.gr, upstream = range[1], downstream = range[2])) ## 400 up and 400 down
  selection <- width(tc.flank) == sum(range) ## only those with width = 800bp
  tc.flank <- tc.flank[selection]
  tc.flank <- tc.flank[order(tc.flank@elementMetadata$interquantile_width)] ## order by IQ width
  # get sequence
  tc.seq <- getSeq(bs_genome, tc.flank)
  #if(remove.na == TRUE){
  #	tc.seq <- clean(tc.seq)
  #}
  return(tc.seq)
}

######### PLOT HEATMAPS
## heatmaps
library(heatmaps)
gen.name <- unname(danRer7CAGEset@sampleLabels)[3]
for(i in 1:length(my.deg.list)){
  tryCatch({
    seqs.tc <- .retrieveSequences(my.deg.list[[i]], range = c(250,250), SeqLengths = seqlengths(BSgenome.Drerio.UCSC.danRer7), bs_genome = BSgenome.Drerio.UCSC.danRer7, remove.na = FALSE)
    s.name <- names(my.deg.list)[i]
    # patterns
    pat <- c("TA","CG", "SS","WW")
    pats <- lapply(pat, function(x){ PatternHeatmap(seqs.tc, x, coords = c(-250,250))})
    patsmooth <- lapply(pats, function(x){ smoothHeatmap(x, output.size=c(1000, 800), sigma = c(5, 3), algorithm="kernel")})
    # save
    png(paste0("~/Desktop/Postdoc/Data_Analysis/CAGEsets/PGC/Heatmaps_Merge_CAGE_RNAseq/Heatmap_250",gen.name,"_",s.name,".png"), height=20, width=40, units="cm", res=150)
    plotHeatmapList(patsmooth, groups=1:4, color=list("Blues","Blues","Blues","Blues"))
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

######### PLOT TBP MOTIF
## plot average TBP motif signal 
data(TBPpwm)

## create a df
operation <- sapply(files, function(x) data.frame(assign(x, matrix(NA, 500))))
df_tbp <- data.frame(lapply(operation, cbind))

for(i in 1:length(my.deg.list)){
  tryCatch({
    seqs.tc <- .retrieveSequences(my.deg.list[[i]], range = c(250,250), SeqLengths = seqlengths(BSgenome.Drerio.UCSC.danRer7), bs_genome = BSgenome.Drerio.UCSC.danRer7, remove.na = FALSE)
    s.name <- names(my.deg.list)[i]
    tbp <- motifScanScores(seqs.tc, motifPWM = TBPpwm)
    tbp <- colMeans(tbp)
    tbp[(length(tbp)+1):500] <- 0
    
    df_tbp[,i] <- tbp
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
# rename columns
for(element in 1:dim(df_tbp)[2]){
  colnames(df_tbp)[element] <- strsplit(colnames(df_tbp)[element], '\\.')[[1]][1]
}


## rescale average signal for plot
library(scales)
library(tidyr)

f <- data.frame(apply(df_tbp, 2, function(x) rescale(x, 0:1)))

df.gg <- f %>% gather(Expression, TBP2Sites)
df.gg$Index <- rep(1:500, times = 8)

ggplot(df.gg, aes(x=Index, y=TBP2Sites, Group=factor(Expression))) +
  geom_line(aes(colour=factor(Expression)), size = 0.1) + ylim(c(35,70)) + 
  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept=0, colour="grey", linetype="longdash") 

######### PLOT META OF WW
dataMetaplot = function(hm_list, binsize=1, colors=gg_col(length(hm_list)), addReferenceLine=FALSE) {
  if (class(hm_list) == "Heatmap") hm_list = list(hm_list) # allow single heatmap argument
  
  if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
    stop("heatmaps must have the same coordinates")
  
  
  if (!length(unique(lapply(hm_list, function(x) xm(x)))) == 1)
    stop("heatmaps must have the same xm values")
  
  n_seq = unique(sapply(hm_list, function(x) x@nseq))
  
  if (!length(n_seq) == 1)
    stop("heatmaps must have the same number of sequences")
  
  coords = hm_list[[1]]@coords
  if (binsize != 1) {
    if (!all(xm(hm_list[[1]]) == 1:width(hm_list[[1]]))) {
      stop("cannot set binsize for heatmaps which are already binned/smoothed")
    }
    breaks = seq(0, width(hm_list[[1]]), by=binsize)
    bin_sums = lapply(hm_list, bin_heatmap, breaks=breaks)
    x_coord = breaks[1:(length(breaks)-1)] + binsize/2 + coords[1]
  } else {
    breaks = xm(hm_list[[1]])
    x_coord = breaks + binsize/2 + coords[1]
    bin_sums = lapply(hm_list, function(x) colSums(image(x)))
  }
  scale_factor = n_seq*(width(hm_list[[1]])/length(x_coord))
  occurrence = lapply(bin_sums, function(x) x/scale_factor)
  max_value = max(vapply(occurrence, max, numeric(1)))
  output = data.frame("x_coord" = x_coord, "occurrence" = occurrence[[1]], "bin_sums" = bin_sums[[1]])
  return(output)
}

#' @importFrom stats aggregate
bin_heatmap = function(hm, breaks) {
  partition = data.frame(pos=xm(hm), value=colSums(image(hm)), bin=cut(xm(hm), breaks))
  aggregate(partition$value, sum, by=list(bin=partition$bin))$x
}

#' @importFrom grDevices hcl
gg_col = function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


## function
na_df <- sapply(files, function(x) data.frame(assign(x, matrix(NA, 150))))
df_ww <- data.frame(lapply(na_df, cbind))

for(i in 1:length(my.deg.list)){
  tryCatch({
    seqs.tc <- .retrieveSequences(my.deg.list[[i]], range = c(150,150), SeqLengths = seqlengths(BSgenome.Drerio.UCSC.danRer7), bs_genome = BSgenome.Drerio.UCSC.danRer7, remove.na = FALSE)
    s.name <- names(my.deg.list)[i]
    # patterns
    pat <- 'TATAA'
    pats <- lapply(pat, function(x){ PatternHeatmap(seqs.tc, x, coords = c(-150,150))})
    meta_l <- dataMetaplot(pats, binsize = 2)
    # fill the dataframe

    df_ww[,i] <- meta_l$occurrence
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

## plot overlap meta
# rename columns
for(element in 1:dim(df_ww)[2]){
  colnames(df_ww)[element] <- strsplit(colnames(df_ww)[element], '\\.')[[1]][1]
}
## keep only interesting groups
df_ww <- df_ww[,c(8,9,12,13)]
df_ww_plot <- df_ww %>% gather(Sample, Occurence)
df_ww_plot$Index <- rep(seq(-149,149, by=2), times = ncol(df_ww))
# plot
ggplot(df_ww_plot, aes(x=Index, y=Occurence, Group=factor(Sample))) +
  geom_line(aes(colour=factor(Sample)), size = 0.6) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept=0, colour="grey", linetype="longdash")


######### PLOT CHROMATIN
## chromatin
setwd('../Data_Analysis/ATACseq/')
files <- list.files(pattern = "\\.danRer7.sorted.goodChr.bam$")
atac.file.names <- files[c(1,2)]
#atac.pgc <- files[grep('PGC', files)]
#atac.soma <- files[grep('somatic', files)][2]
## load transcipt file
txs <- loadDb('../annotation/txdb_DanRer7.sqlite')
blacklist <- read.csv('../danRer10_blacklist.bed', sep = '\t', header = F)
blacklist <- GRanges(seqnames = blacklist$V1, IRanges(start = blacklist$V2, end = blacklist$V3))

call.atac <- lapply(atac.file.names, function(x){
  atac.file <- readGAlignments(x)
  atacR = granges(atac.file, use.names=T, use.mcols=F)
  cut0 = GRanges(seqnames(atacR), IRanges(start(atacR)+5, start(atacR)+5), strand(atacR))
  cut1 = GRanges(seqnames(atacR), IRanges(end(atacR)-4, end(atacR)-4), strand(atacR))
  a.pgc = GRanges(seqnames(atacR), IRanges(start(cut0), end(cut1)), strand(atacR)) 
})

## calculate fold change obs coverage vs expected
cov.function <- function(a, resize = 1){
  #a <- subset(a, width(a) < 120)
  seqlevels(a) = seqlevels(BSgenome.Drerio.UCSC.danRer7)
  seqinfo(a)   = seqinfo(BSgenome.Drerio.UCSC.danRer7)
  a = trim(a)
  ## remove blacklisted regions
  a <- subsetByOverlaps(a, blacklist, invert = T)
  ## Tn5 cut sites: resize = 1
  a <- resize(a, width = resize, fix = "start", ignore.strand = FALSE)
  
  a_plus <- a[strand(a) == '+']
  a_minus <- a[strand(a) == '-']
  
  genome.size <- 1.42e9
  expected.cov <- sum(width(a))/genome.size
  ## select only nucleosomes
  FoldChange.minus <- GRanges(coverage(a_minus) / expected.cov, strand = '-')
  FoldChange.plus <- GRanges(coverage(a_plus) / expected.cov, strand = '+')
  FoldChange.strand <- c(FoldChange.minus, FoldChange.plus)
  
  #export.bed(FoldChange, paste0(x, '.bed'))
  return(FoldChange.strand)
}
atac.fold.change.pgc <- cov.function(call.atac[[1]], resize = 20)
atac.fold.change.soma <- cov.function(call.atac[[2]], resize = 20)

# plot
library(genomation)

df_chromatin <- sapply(files, function(x) data.frame(assign(x, matrix(NA, 500))))
df_chromatin <- data.frame(lapply(df_chromatin, cbind))
  
for(i in 1:length(my.deg.list)){
  tryCatch({
    f <- GRanges(seqnames = seqnames(my.deg.list[[i]]), IRanges(ranges(my.deg.list[[i]])),
                 strand = strand(my.deg.list[[i]]))
    f <- resize(f, 500, fix = 'center')
    sm <- ScoreMatrix(target = atac.fold.change.soma, windows = f, weight.col = 'score', 
                      strand.aware = T)
    col_data <- data.frame(colMeans(sm))
    
    df_chromatin[,i] <- col_data
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

# rename columns
for(element in 1:dim(df_chromatin)[2]){
  colnames(df_chromatin)[element] <- strsplit(colnames(df_chromatin)[element], '\\.')[[1]][1]
}

## rescale average signal for plot
library(scales)
library(tidyr)

df_chromatin_plot <- data.frame(apply(df_chromatin, 2, function(x) rescale(x, 0:1)))
df_chromatin_plot <- df_chromatin_plot[,c(8,9,12,13)]
df.gg <- df_chromatin_plot %>% gather(Expression, Tn5CutSites)
df.gg$Index <- rep(seq(-249,250, by = 1), times = ncol(df.gg))

ggplot(df.gg, aes(x=Index, y=Tn5CutSites, Group=factor(Expression))) +
  geom_line(aes(colour=factor(Expression)), size = 0.6) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept=0, colour="grey", linetype="longdash")


