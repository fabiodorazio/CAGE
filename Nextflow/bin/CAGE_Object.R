library(CAGEr)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)

# inputs mart and flanks around Ensembl TSS to create promoters from
# promoter data.frame is returned, duplicates based on same start and end are removed
promoters_from_cage <- function (CAGEset) {

get_CAGE_TCs <- function (CAGEset,
                          sample,
                          min_tpm = 1,
                          librarySizes = library_sizes, 
                          clustering_method = "powerLaw",
                          fitInRange = c(5,1000),
                          qLow = 0.1,
                          qUp = 0.9,
                          multicore = FALSE)

  
  reference.library.size <- 10^floor(log10(median(librarySizes))) # modified by Anja (fails on tag sum > max. integer)
  reference.alpha <- -1*median(apply(CAGEset@tagCountMatrix, 2, function(x) {
    CAGEr:::.fit.power.law.to.reverse.cumulative(values = as.integer(x), val.range = fitInRange)
  })[1,])
  message("Alpha for ",sample, ": ", reference.alpha)
  # save reverse cumulative plot
  pdf("../Outputs/mm9_pooledCTSS_revcumulative.pdf")
  plotReverseCumulatives(CAGEset, fitInRange = fitInRange, onePlot = TRUE)
  dev.off()
  
  normalizeTagCount(CAGEset,
                    method = "powerLaw",
                    fitInRange = fitInRange,
                    alpha = reference.alpha,
                    T = reference.library.size)
  
  ## cluster cage tags in CTSS
  clusterCTSS(CAGEset,
              threshold = min_tpm,
              thresholdIsTpm = TRUE,
              method = "distclu",
              useMulticore = as.logical(multicore), nrCores = multicore)
  ## determines promoter width 
  cumulativeCTSSdistribution(CAGEset, clusters = "tagClusters",
                             useMulticore = as.logical(multicore), nrCores = multicore)
  quantilePositions(CAGEset, clusters = "tagClusters", qLow = qLow, qUp = qUp,
                    useMulticore = as.logical(multicore), nrCores = multicore)
  tc_df <- tagClusters(CAGEset, sample = sample,
                       returnInterquantileWidth = TRUE, qLow = qLow, qUp = qUp)
  return(tc_df)
}


### run functions

promoters_df <- promoters_from_cage(CAGEset = args[1])
