library(DESeq2)
library(biomaRt)
library(ChIPpeakAnno)
library(ChIPseeker)
library(AnnotationDbi)
library(BSgenome.Drerio.UCSC.danRer7)
library(seqPattern)
library(CAGEr)
source("~/Desktop/Postdoc/R_scripts/gcCAGE/R/CAGE_Functions.R")
source("../src/Run_DEseq2.R")

## plot the open chromatin on germ cell genes
## get open chromatin signal
## generate big data frame where germ cell genes are separated from somatic genes and ubiquitous

args<-commandArgs(TRUE)
file_path <- args[1]
#setwd('perGeneCounts_Stages/')
setwd(file_path)

files = list.files(pattern = "\\.tab")

readTagsPerGene <- function(x){
  out <- read.csv(x, skip = 4, header = FALSE, sep = "\t", row.names = 1,
                  col.names = c("","totalCount", "forwardCount", "reverseCount"))
  out
}

allReadCount <- lapply(files, readTagsPerGene)
totalReadCount <- sapply(allReadCount, function(x) x$totalCount)
rownames(totalReadCount) <- rownames(allReadCount[[1]])
## add column names
spl_names <- strsplit(files, "\\_R")
s <- as.character(lapply(spl_names, function(x) x[[1]][1]))
colnames(totalReadCount) <- s

## differentially expressed genes at prim5

## subset stages
prim5 <- totalReadCount[,c(17:20)]
condition <- rep(c('PGC', 'Soma'), times = 2)

dds <- run.deseq.genes(prim5, condition)
res <- results(dds)
res_final <- data.frame(res[,c(2,6)])
# load TBP2 res object
res.tbp2 <- read.csv("../data/Res_TBP2_newMM.txt", sep = '\t')

prepare_res <- function(x){
  res.data <- data.frame(x)
  res.data$GeneType <- 'Ubiquitous'
  res.data <- na.omit(res.data)
  
  res.data$GeneType[(res.data$padj < 0.05) & (res.data$log2FoldChange > 1)] <- 'Soma'
  res.data$GeneType[(res.data$padj < 0.05) & (res.data$log2FoldChange < -1)] <- 'PGC'
  res.data <- data.frame(row.names = rownames(res.data), "GeneType"=res.data$GeneType)
  return(res.data)
}
res.data.wt <- prepare_res(res)
res.data.tbp2 <- prepare_res(res.tbp2)

## get expression data

tpm <- read.csv('../tpmStagesFiltered.txt', sep = '\t')

tpm$TpmMeanPGChigh <- rowMeans(tpm[c("sHighPGC1", "sHighPGC2")],)
tpm$TpmMeanSomahigh <- rowMeans(tpm[c("sHighSoma1", "sHighSoma2")],)
tpm$TpmMeanPGCprim <- rowMeans(tpm[c("sPrim5PGC1", "sPrim5PGC2")],)
tpm$TpmMeanSomaprim <- rowMeans(tpm[c("sPrim5Soma1", "sPrim5Soma2")],)
# subset the mean of replicates
tpm_final <- tpm[,c((ncol(tpm)-3):(ncol(tpm)))]
# merge with Diff expr analysis
final_data_PGC_1 <- merge(tpm_PGC, res.data.wt, by = 0)
final_data_TBP_1 <- merge(tpm_PGC, res.data.tbp2, by = 0)

## calculate the Expression variation
tpm_final$Expr_var_PGC <- tpm_final$TpmMeanPGCprim/tpm_final$TpmMeanPGChigh
tpm_final$Expr_var_Soma <- tpm_final$TpmMeanSomaprim/tpm_final$TpmMeanSomahigh

##annotate
mart <- useMart(host='dec2017.archive.ensembl.org',
                biomart='ENSEMBL_MART_ENSEMBL',
                dataset = "drerio_gene_ensembl")
#listAttributes(mart)
fish_genes <- getBM(attributes = c("ensembl_gene_id",
                                   "transcript_count",
                                   "transcript_length",
                                   "strand",
                                   "gene_biotype"),
                    mart = mart)

## merge with RNAseq
final_data_PGC_2 <- merge(final_data_PGC_1, fish_genes, by.x = "Row.names", by.y = "ensembl_gene_id")
final_data_TBP_2 <- merge(final_data_TBP_1, fish_genes, by.x = "Row.names", by.y = "ensembl_gene_id")

###########
## merge with cage to determine promoter class (sharp, broad, no prom)

myCAGEset <- readRDS("../data/CAGEset_PGC_soma_Early_Late.rds")

## shift/no shift
aggregateTagClusters(myCAGEset, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters")
quantilePositions(myCAGEset, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9, useMulticore = FALSE)
## shifting promoters
myCAGEset_PGC <- myCAGEset
myCAGEset_Soma <- myCAGEset

scoreShift(myCAGEset_PGC, groupX = "PGC_high_r2_1_S3", groupY = "PGC_prim5_r4_1_S5",
           testKS = TRUE, useTpmKS = FALSE)
scoreShift(myCAGEset_Soma, groupX = "soma_high_r2_1_S4", groupY = "soma_prim5_r4_1_S6",
           testKS = TRUE, useTpmKS = FALSE)

shifting.promoters_PGC <- getShiftingPromoters(myCAGEset_PGC,
                                               tpmThreshold = 5, scoreThreshold = 0.6,
                                               fdrThreshold = 0.01)
shifting.promoters_PGC$Shifting <- "PGC_shift"


## annotate shifting promoters

txdb <- loadDb("../annotation/txdb_DanRer7.sqlite")
shifting.promoters_PGC_anno <- annotate.cage.peaks(shifting.promoters_PGC)
shifting.promoters_PGC_anno <- shifting.promoters_PGC_anno[,c("geneId","Shifting")]

## define promoter width per gene
cage_names <- unname(myCAGEset@sampleLabels[c(1,3,4,6)])

tc_list <- lapply(cage_names, function(x){
  tc <- tagClusters(myCAGEset, sample = x, returnInterquantileWidth = TRUE,
                    qLow = 0.1, qUp = 0.9)
  tc <- tc[tc$chr %in% paste("chr",1:25, sep = ""),]
  tc$sampleID = x
  tc$Width_Class <- "Broad"
  tc$Width_Class[tc$interquantile_width <= 9] <- "Sharp"
  return(tc)
})
names(tc_list) <- cage_names

# get position of TATA boxes centered to PGC or somatic TSSs at prim5
data(TBPpwm)

.tata_positions <- function(cage_obj){
  cage_obj_range <- GRanges(seqnames = cage_obj$chr, IRanges(start = cage_obj$dominant_ctss,
                                                             end = cage_obj$dominant_ctss),
                            strand = cage_obj$strand)
  values(cage_obj_range) <- cage_obj[,-c(1:5)]
  cage_obj_range <- resize(cage_obj_range, 150, fix = 'end')
  # get sequences
  cage_seqs <- getSeq(Drerio, cage_obj_range)
  
  m <- motifScanHits(regionsSeq =cage_seqs, motifPWM = TBPpwm, minScore = "60%")
  
  max_duplicates <- ddply(m,.(sequence), nrow)
  max_duplicates <- subset(max_duplicates, max_duplicates$V1 <= 6)
  max_duplicates <- max(max_duplicates$V1)
  
  # create columns to populate
  num_col_original <- ncol(m) + 1
  new_names <- paste0('TATA_n_',2:max_duplicates)
  m[,c(as.numeric(num_col_original):as.numeric(ncol(m)+max_duplicates-1))] <- 0
  colnames(m)[num_col_original:ncol(m)] <- new_names
  # fill the df with 6 empty rows to avoid data loss
  fake_df <- head(m)
  fake_df[,1] <- "NA"
  fake_df[,2:ncol(fake_df)] <- 0
  m <- rbind(fake_df, m)
  
  
  for(element in 2:nrow(m)){
    tryCatch({
      if(m$sequence[element] == m$sequence[element - 1] & m$sequence[element] != m$sequence[element - 2]){
        m$TATA_n_2[element-1] <- m$position[element]
      }
      else if(m$sequence[element] == m$sequence[element - 2] & m$sequence[element] != m$sequence[element - 3]){
        m$TATA_n_3[element-2] <- m$position[element]
      }
      else if(m$sequence[element] == m$sequence[element - 3] & m$sequence[element] != m$sequence[element - 4]){
        m$TATA_n_4[element-3] <- m$position[element]
      }
      else if(m$sequence[element] == m$sequence[element - 4] & m$sequence[element] != m$sequence[element - 5]){
        m$TATA_n_5[element-4] <- m$position[element]
      }
      else if(m$sequence[element] == m$sequence[element - 5] & m$sequence[element] != m$sequence[element - 6]){
        m$TATA_n_6[element-5] <- m$position[element]
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  # remove top rows
  m <- m[-c(1:6),]
  # count repeated TATAs
  m3 <- m %>% 
    group_by(as.numeric(sequence)) %>%
    dplyr::summarise(N = n()) %>%
    data.frame()
  
  # remove duplicated rows
  n <- m[!duplicated(m$sequence),]
  
  # add count to original table
  n$num_TATA <- m3$N
  # reconsitute original table
  n$value <- NULL
  colnames(n)[2] <- 'TATA_n_1'
  
  # average distance
  n1 <- n
  n1[n1 == 0] <- NA
  n1$Avg_TATA_dist = apply(n1,1,function(x) mean(as.numeric(x[2:7])[!is.na(x[2:7])]))
  
  # add to original
  n$Avg_TATA_dist <- n1$Avg_TATA_dist
  
  # merge
  cage_obj_tata <- merge(cage_obj, n, by.x = "cluster", by.y = "sequence", all = TRUE)
  
}

# apply to tc_list
tc_list1 <- lapply(tc_list, .tata_positions)

# add ENSEMBL IDs
tc_list_anno <- lapply(tc_list1, annotate.cage.peaks)
tc_list_anno1 <- lapply(tc_list_anno, function(x){
  tag_cluster <- x[,-c(33:35)]
  return(tag_cluster)
})
