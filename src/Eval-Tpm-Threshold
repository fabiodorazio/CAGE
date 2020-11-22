library(BSgenome.Drerio.UCSC.danRer10)
library(CAGEr)
library(GenomicRanges)

paths <- grep(list.files(), pattern = "bowtiemapping.bam", value = T)

names <- c("high","prim5")
paths <- paths[c(6,12)]

#####################################################
### Run CAGEr pipeline to get consensus clusters ####

danRer10CAGEset <- new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer10", inputFiles = paths, inputFilesType = "bam", 
                       sampleLabels = names)
getCTSS(danRer10CAGEset)

## check correlation between samples
run_cage_r <- function(myCAGEset2){
  normalizeTagCount(myCAGEset2, method = "powerLaw",
                    fitInRange = c(5, 1000), alpha = 1.18, T = 1*10^6)

  
  clusterCTSS(object = myCAGEset2, threshold = 1, thresholdIsTpm = TRUE,
              nrPassThreshold = 1, method = "distclu", maxDist = 20,
              removeSingletons = TRUE, keepSingletonsAbove = 5)
  cumulativeCTSSdistribution(myCAGEset2, clusters = "tagClusters")
  quantilePositions(myCAGEset2, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
  aggregateTagClusters(myCAGEset2, tpmThreshold = 5,
                       qLow = 0.1, qUp = 0.9, maxDist = 100)
  cumulativeCTSSdistribution(myCAGEset2, clusters = "consensusClusters")
  return(myCAGEset2)
}

CAGE_stage <- run_cage_r(danRer10CAGEset)
CAGE_PGC <- run_cage_r(myCAGEsetPGC)

list_cage  <- c(CAGE_stage, CAGE_PGC)


evaluate_threshold_quality <- function(x){
  list1 <- list()
  for(cage_obj in seq_along(x)){
    xq <- x[[cage_obj]] ## cage object
    names <- unname(xq@sampleLabels)
    for(n in names){
      ## calculate the ratio between n of cage tags and tag clusters
      n_cage_tags <- sum(xq@tagCountMatrix[n])
      n_tag_clusters <- nrow(xq@tagClusters[[n]])
      eval_ratio <- n_cage_tags/n_tag_clusters
      
      list1[[n]] <- eval_ratio
    }

  }  
  return(list1)
}
