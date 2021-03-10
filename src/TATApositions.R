# requires:
# - cager object
# - TSS coordinates

# get position of TATA boxes
library(BSgenome.Drerio.UCSC.danRer7)
library(seqPattern)
data(TBPpwm)

cage_obj_range <- GRanges(seqnames = cage_obj$chr, IRanges(start = cage_obj$dominant_ctss,
                                                           end = cage_obj$dominant_ctss),
                          strand = cage_obj$strand)
values(cage_obj_range) <- cage_obj[,-c(1:5)]
cage_obj_range <- resize(cage_obj_range, 150, fix = 'end')
# get sequences
cage_seqs <- getSeq(Drerio, cage_obj_range)

m <- motifScanHits(regionsSeq =cage_seqs, motifPWM = TBPpwm, minScore = "90%")

library(plyr) # to count max duplicates
max_duplicates <- ddply(m,.(sequence), nrow)
max_duplicates <- subset(max_duplicates, max_duplicates$V1 <= 6)
max_duplicates <- max(max_duplicates$V1)

# create columns to populate
num_col_original <- ncol(m) + 1
new_names <- paste0('TATA_n_',2:max_duplicates)
m[,c(as.numeric(num_col_original):as.numeric(ncol(m)+max_duplicates-1))] <- 0
colnames(m)[num_col_original:ncol(m)] <- new_names


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

m2 <- m
m2$test <- paste0('t', m2$sequence)
# count repeated TATAs
m3 <- m %>% 
  group_by(sequence) %>% 
  dplyr::summarise(N = n()) %>%
  data.frame()

# remove duplicated rows
n <- m[!duplicated(m$sequence),]
# add count to original table
n$num_TATA <- m3$N

# reconsitute original table
n$value <- NULL
colnames(n)[2] <- 'TATA_n_1'
# merge
cage_obj_tata <- merge(cage_obj, n, by.x = "cluster", by.y = "sequence", all = TRUE)
