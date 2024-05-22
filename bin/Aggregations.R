## merge with shifting promoters
# 1) shifting.promoters_PGC_anno
## load and retrieve shifting scores
library(AnnotationDbi)
library(ChIPpeakAnno)
library(ChIPseeker)
library(reshape2)
library(ggplot2)
source("../src/CAGE_Functions.R")

output_path <- '../Output'
setwd('../data')

# load datasets
myCAGEset_PGC <- readRDS('CAGEsetPGC_Shift_HighToPrim5.rds')
myCAGEset_Soma <- readRDS('CAGEsetSoma_Shift_HighToPrim5.rds')

shifting.promoters_PGC <- getShiftingPromoters(myCAGEset_PGC,
                                               tpmThreshold = 5, scoreThreshold = 0.6,
                                               fdrThreshold = 0.01)
shifting.promoters_PGC$Shifting <- "PGC_shift"

shifting.promoters_Soma <- getShiftingPromoters(myCAGEset_Soma,
                                                tpmThreshold = 5, scoreThreshold = 0.6,
                                                fdrThreshold = 0.01)
shifting.promoters_Soma$Shifting <- "Soma_shift"
# add ensembls
txdb <- loadDb("annotation/txdb_DanRer7.sqlite")

shifting.promoters_PGC_anno <- annotate.cage.peaks(shifting.promoters_PGC)
shifting.promoters_PGC_anno <- shifting.promoters_PGC_anno[,c("geneId","Shifting", "shifting.score")]

shifting.promoters_Soma_anno <- annotate.cage.peaks(shifting.promoters_Soma)
shifting.promoters_Soma_anno <- shifting.promoters_Soma_anno[,c("geneId","Shifting", "shifting.score")]

#### create table to merge
s_class <- merge(shifting.promoters_PGC_anno, shifting.promoters_Soma_anno, by = 'geneId', all = TRUE)
# remove NA
s_class$Shifting.x[is.na(s_class$Shifting.x)] <- 'No_shift'
s_class$Shifting.y[is.na(s_class$Shifting.y)] <- 'No_shift'
s_class[is.na(s_class)] <- 0

s_class$Shift_Prim5[(s_class$Shifting.x == 'PGC_shift') & (s_class$Shifting.y == 'No_shift')] <- "PGC"
s_class$Shift_Prim5[(s_class$Shifting.x == 'PGC_shift') & (s_class$Shifting.y == "Soma_shift")] <- "Both"
s_class$Shift_Prim5[(s_class$Shifting.x == 'No_shift') & (s_class$Shifting.y == "No_shift")] <- "None"
s_class$Shift_Prim5[(s_class$Shifting.x == 'No_shift') & (s_class$Shifting.y == "Soma_shift")] <- "Soma"

# final table to merge
s <- s_class[,c(1,3,5,6)]
# rename columns
colnames(s)[2] <- 'Shifting_Score_prim5_PGC'
colnames(s)[3] <- 'Shifting_Score_prim5_Soma'

####
merge_3 <- function(x){
  # add shifting class
  my_m <- merge(x, s, by = 'geneId', all = TRUE)
  #my_m$Shifting_Score_prim5_PGC[is.na(my_m$Shifting_Score_prim5_PGC)] <- 0
  #my_m$Shifting_Score_prim5_Soma[is.na(my_m$Shifting_Score_prim5_Soma)] <- 0
  my_m$Shift_Prim5[is.na(my_m$Shift_Prim5)] <- 'None'
  
  return(my_m)
}
##### 3
final_data_3_list <- lapply(tc_list_anno1, merge_3)

## merge CAGE with RNAseq
merge_4 <- function(x){
  my_m <- merge(final_data_PGC_2, x, by.x = "Row.names", by.y = "geneId",
        all.x = TRUE)
  # add shifting class
  my_m <- merge(my_m, s_class[,c(1,4)], by.x = 'Row.names', by.y = 'geneId', all = TRUE)
  return(my_m)
}
##### 4
final_data_PGC_4 <- lapply(final_data_3_list, merge_4)
######
final_data_TBP_4 <- merge(final_data_TBP_2, final_data_PGC_3, by.x = "Row.names", by.y = "geneId",
                          all.x = TRUE)

######## 5
## merge with atac seq to count elements (count enhancers)
setwd('../Summit_MACS/')
file_summit <- list.files()[grep('summit', list.files())]

.import.summit <- function(f_summit){
  summit <- read.csv(f_summit, sep = '\t', header = F)[,c(1:3)]
  colnames(summit) <- c('seqnames', 'start', 'end')
  s_range <- GRanges(seqnames = summit$seqnames, IRanges(start = summit$start, end = summit$start))
  return(s_range)
}
pgc.summit <- .import.summit(file_summit[1])
soma.summit <- .import.summit (file_summit[2])

## annotate peaks to the nearest gene
ensembl.ids.atac <- function(grange.obj){
  anno.ensembl.id <- annotatePeakInBatch(grange.obj, AnnotationData=toGRanges(txdb), 
                                         output="nearestBiDirectionalPromoters",
                                         bindingRegion=c(-10000, 10000))
  anno.ensembl.id <- data.frame("geneId" = anno.ensembl.id$feature, "distanceFromProm" = anno.ensembl.id$distance,
                                "insideFeature" = anno.ensembl.id$insideFeature)
  anno.ensembl.id <- anno.ensembl.id[order(anno.ensembl.id$geneId),]
  return(anno.ensembl.id)
  
}

all.peaks.pgc <- ensembl.ids.atac(pgc.summit)
all.peaks.soma <- ensembl.ids.atac(soma.summit)

atac_list <- list(all.peaks.pgc, all.peaks.soma)

## get peaks coordinates for each promoter
.atac_positions <- function(y){
  max_duplicates <- ddply(y,.(geneId), nrow)
  max_duplicates <- max(max_duplicates$V1)
  
  # create columns to populate
  num_col_original <- ncol(y) + 1
  new_names <- paste0('Peak_n_',2:max_duplicates)
  y[,c(as.numeric(num_col_original):as.numeric(ncol(y)+max_duplicates-1))] <- 0
  colnames(y)[num_col_original:ncol(y)] <- new_names
  # ascending order
  y <- y[order(y$geneId),]
  # sort in ascending order within groups
  y <- y %>%
    group_by(geneId) %>%
    arrange(geneId, distanceFromProm) %>%
    data.frame()
  
  # fill the df with 6 empty rows to avoid data loss
  fake_df <- head(y)
  fake_df[,1] <- "f_ENSEMBL"
  fake_df[,2:ncol(fake_df)] <- 0
  y <- rbind(fake_df, y)
  
  # fill the data frame
  for(element in 1:nrow(y)){
    tryCatch({
      if(y$geneId[element] == y$geneId[element - 1] & y$geneId[element] != y$geneId[element - 2]){
        y$Peak_n_2[element-1] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 2] & y$geneId[element] != y$geneId[element - 3]){
        y$Peak_n_3[element-2] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 3] & y$geneId[element] != y$geneId[element - 4]){
        y$Peak_n_4[element-3] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 4] & y$geneId[element] != y$geneId[element - 5]){
        y$Peak_n_4[element-4] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 5] & y$geneId[element] != y$geneId[element - 6]){
        y$Peak_n_4[element-5] <- y$distanceFromProm[element]
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  # remove top rows
  y <- y[-c(1:6),]
  # count repeated Peaks
  y3 <- y %>% 
    group_by(geneId) %>%
    dplyr::summarise(N = n()) %>%
    data.frame()
  
  # remove duplicated rows
  z <- y[!duplicated(y$geneId),]
  
  # add count to original table
  z$num_ATAC_peaks <- y3$N
  # reconsitute original table
  z$insideFeature <- NULL
  colnames(z)[2] <- 'Peak_n_1'
  
  # average distance
  z1 <- z
  z1[z1 == 0] <- NA
  z1$Avg_Peak_dist = apply(z1,1,function(x) mean(as.numeric(x[2:7])[!is.na(x[2:7])]))
  
  # add to original
  z$Avg_Peak_dist <- z1$Avg_Peak_dist
  
  return(z)
  
  #cage_obj_tata <- merge(cage_obj, z, by.x = "cluster", by.y = "sequence", all = TRUE)
}

atac_list_positions <- lapply(atac_list, .atac_positions)

## merge atac info
final_data_PGC_5 <- merge(final_data_3_list[[1]], atac_list_positions[[1]], by = 'geneId', all.x = TRUE)
final_data_Soma_5 <- merge(final_data_3_list[[2]], atac_list_positions[[2]], by = 'geneId', all = TRUE)

# save
final_data_5_list <- list(final_data_PGC_5,final_data_Soma_5)
names(final_data_5_list) <- c('PGC_prim5', 'Soma_prim5')

## get mean of distance from promoter for same ENSEMBL IDs
library(tidyverse)
all.peaks.pgc.comp <- all.peaks.pgc %>% 
  group_by(geneId) %>% 
  summarise(Distance_ATAC_Mean = mean(distanceFromProm), Num_ATAC_peaks = n()) %>%
  data.frame()
all.peaks.soma.comp <- all.peaks.soma %>% 
  group_by(geneId) %>% 
  summarise(Distance_ATAC_Mean = mean(distanceFromProm), Num_ATAC_peaks = n()) %>%
  data.frame()

final_data_PGC_5 <- merge(final_data_PGC_4, all.peaks.pgc.comp, by.x = "Row.names", by.y = "geneId")
final_data_TBP_5 <- merge(final_data_TBP_4, all.peaks.pgc.comp, by.x = "Row.names", by.y = "geneId")

## CpG islands
library(AnnotationDbi)
library(ChIPpeakAnno)
library(rtracklayer)
source("../src/CAGE_Functions.R")

cpg_island <- import.bed('CpG.prediction.Christopher.UCSC.track.bed')
# annotate Cpg
anno.ensembl.id <- annotatePeakInBatch(cpg_island, AnnotationData=toGRanges(txdb), 
                                       output="nearestBiDirectionalPromoters",
                                       bindingRegion=c(-1000, 1000))
anno.ensembl.id <- data.frame("geneId" = anno.ensembl.id$feature, "distanceFromProm" = anno.ensembl.id$distance)

library(plyr)
library(dplyr)
.cpg_positions <- function(y){
  max_duplicates <- ddply(y,.(geneId), nrow)
  max_duplicates <- max(max_duplicates$V1)
  
  # create columns to populate
  num_col_original <- ncol(y) + 1
  new_names <- paste0('Cpg_n_',2:max_duplicates)
  y[,c(as.numeric(num_col_original):as.numeric(ncol(y)+max_duplicates-1))] <- 0
  colnames(y)[num_col_original:ncol(y)] <- new_names
  # ascending order
  y <- y[order(y$geneId),]
  # sort in ascending order within groups
  y <- y %>%
    group_by(geneId) %>%
    arrange(geneId, distanceFromProm) %>%
    data.frame()
  
  # fill the df with 6 empty rows to avoid data loss
  fake_df <- head(y)
  fake_df[,1] <- "f_ENSEMBL"
  fake_df[,2:ncol(fake_df)] <- 0
  y <- rbind(fake_df, y)
  
  # fill the data frame
  for(element in 1:nrow(y)){
    tryCatch({
      if(y$geneId[element] == y$geneId[element - 1] & y$geneId[element] != y$geneId[element - 2]){
        y$Cpg_n_2[element-1] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 2] & y$geneId[element] != y$geneId[element - 3]){
        y$Cpg_n_3[element-2] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 3] & y$geneId[element] != y$geneId[element - 4]){
        y$Cpg_n_4[element-3] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 4] & y$geneId[element] != y$geneId[element - 5]){
        y$Cpg_n_5[element-4] <- y$distanceFromProm[element]
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  # remove top rows
  y <- y[-c(1:6),]
  # count repeated TATAs
  y3 <- y %>% 
    group_by(geneId) %>%
    dplyr::summarise(N = n()) %>%
    data.frame()
  
  # remove duplicated rows
  n <- y[!duplicated(y$geneId),]
  
  # add count to original table
  n$num_CpG <- y3$N
  # reconsitute original table
  colnames(n)[2] <- 'Cpg_n_1'
  
  return(n)

}

cpg_all <- .cpg_positions(anno.ensembl.id)
  
final_data_3_list_Cpg <- lapply(final_data_3_list, function(x) merge(x, cpg_all, by = 'geneId'))

############# THIS IS THE FINAL TABLE #################

.clean_final_table <- function(x){
  f_5 <- final_data_5_list[[1]]
  f_5 <- f_5[,c(1,(ncol(f_5)-7):ncol(f_5))]
  # remove no shifting
  sub <- subset(x, x$Shift_Prim5 != 'None')
  # merge with ATAC PGC
  sub <- merge(sub, f_5, by = 'geneId')
  # merge with expression data
  sub_expr <- merge(sub, tpm_final, by.x = 'geneId', by.y = 0)
}
final_data_6_list <- lapply(final_data_3_list_Cpg, .clean_final_table)

.merge_atac <- function(x){
  x <- x[!duplicated(x[,c('geneId')]),]
  my_atac <- merge(x, atac_list_positions[[1]], by = 'geneId')
  my_atac <- merge(my_atac, atac_list_positions[[2]], by = 'geneId')
}
final_data_3_list_Atac <- lapply(final_data_3_list_no_Atac, .merge_atac)

# GO
library(clusterProfiler)
library(org.Dr.eg.db)
f <- subset(final_data_3_list[[2]], final_data_3_list[[2]]$Shift_Prim5 == 'PGC')
ego <- enrichGO(f$geneId, 
                     OrgDb = org.Dr.eg.db, keyType = 'ENSEMBL', ont = 'BP', pvalueCutoff = 0.1, 
                     universe = rownames(tpm))
barplot(ego_down, showCategory = 25, font.size = 8)

# check outliers
setwd(paste0(output_path, '/Plots/'))
 # binarize & fill
.binarize <- function(b){
  # fill missing values in Shift_Class
  b$Shift_Class[is.na(b$Shift_Class)] <- 'None'
  # binarize Width_Class
  b$Width_Class <- ifelse(b$Width_Class == 'Broad', 1, 0)
  b$Shifting_Prom <- ifelse(b$Shifting_Prom == 'No_Shift', 1, 0)
  return(b)
}
final_data_5_list <- lapply(final_data_5_list, .binarize)
# relevant columns
dd <- final_data_5_list[[1]][,c(5,8,10,11,12,13,14,16:24,29,32,32:41)]

# plot boxplots
pdf('boxplots_Width_Class.pdf', width=10, height=20)
ggplot(melt(dd, id.vars = "Shifting_Prom"),
       aes(x = variable, y = value, col = Shifting_Prom)) +
  geom_violin() + facet_wrap( ~variable, ncol = 2, scales = "free") +
  theme_bw()
dev.off()
# plot densities
pdf('densities.pdf', width=10, height=40)
ggplot(melt(dd, id.vars = "Shifting_Prom"),
       aes(x = value, col = Shifting_Prom)) +
  geom_density() + facet_wrap( ~variable, ncol = 1, scales = "free") +
  theme(legend.position = "bottom") + theme_bw()
dev.off()
#...............

## prepare the data for ML
f <- function(x){
  # convert NA to zeros
  #x[is.na(x$Shifting_Prom)] <- NULL
  # remove purposeless columns and arrange
  #x <- x[,c(35,5,8,10,11,12,13,14,16:24,29,32:34,36:42)]
  x <- x[,-which(names(x) %in% c('geneId', 'seqnames', 'start', 'end', 'strand', 'cluster', 'dominant_ctss',
             'q_0.1', 'q_0.9', 'sampleID', 'annotation', 'geneChr', 'geneStart', 'geneEnd',
             'geneStrand', 'transcriptId', 'Shift_Prim5'))]
  # binarize
  x$Width_Class <- ifelse(x$Width_Class == 'Broad', 1, 0)
  # convert NA to zeros
  x[is.na(x)] <- 0
  # remove outliers 
  outliers <- function(x) {
    
    Q1 <- quantile(x, probs=.25)
    Q3 <- quantile(x, probs=.75)
    iqr = Q3-Q1
    
    upper_limit = Q3 + (iqr*1.5)
    lower_limit = Q1 - (iqr*1.5)
    
    x > upper_limit | x < lower_limit
  }
  remove_outliers <- function(df, cols = names(df)) {
    for (col in cols) {
      df <- df[!outliers(df[[col]]),]
    }
    df
  }
  x <- remove_outliers(x, colnames(x)[c(2:5,8,17)])
}

r <- lapply(final_data_5_list, f)

## plot correlations
r1 <- r[[2]]
r1$Width_Class <- NULL

library(tidyverse)

for(n in 1:length(r)){
  my_names <- names(r)[n]
  
  r.gathered <- r[[n]] %>%
    as_data_frame() %>%
    gather(key = "variable", value = "value",
           -shifting.score, -Shifting_Prom)
  pdf(paste0(output_path, '/Corr_', my_names, '.pdf'))
  ggplot(r.gathered, aes(x = value, y = shifting.score)) +
    geom_point(aes(color = Shifting_Prom)) +
    facet_wrap(~variable)+
    scale_color_viridis_d()
  dev.off()
}


# boxplot
dd2 <- r[[1]]

# plot boxplots
pdf('Shift_Classes_boxplots_no_outliers.pdf', width=10, height=20)
ggplot(melt(dd2, id.vars = "Shift_Class"),
       aes(x = variable, y = value, col = Shift_Class)) +
  geom_violin() + facet_wrap( ~variable, ncol = 2, scales = "free") +
  theme_bw()
dev.off()

# plot densities
pdf('densities_no_outliers.pdf', width=10, height=40)
ggplot(melt(dd2, id.vars = "Shifting_Prom"),
       aes(x = value, col = Shifting_Prom)) +
  geom_density() + facet_wrap( ~variable, ncol = 1, scales = "free") +
  theme(legend.position = "bottom") + theme_bw()
dev.off()

# scaling
library("ggplot2")
library("reshape2")
library("dplyr")
library("MASS")

set.seed(123)
dd2 <- dd2[dd2$Shift_Class != 'None',]
#dd2 <- dd2[,-c(9,18)]

training.samples <- dd2$Shift_Class %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- dd2[training.samples, ]
test.data <- dd2[-training.samples, ]

preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

model <- lda(Shift_Class~., data = train.transformed)
lda.data <- cbind(train.transformed, predict(model)$x)
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = Shift_Class))

# make predictions
predictions <- model %>% predict(test.transformed)
# how well the model predicts the classes
mean(predictions$class==test.transformed$Shift_Class)


#### remove non meaningfull info
dd3 <- r[[1]]
dd2 <- dd3[dd3$Shift_Class != 'None',]
dd2 <- dd2[,c(1:19)]
dd2 <- dd2[,c(1,2,4,5,10,17)]
# plot boxplots
pdf('boxplots_no_outliers_scaled.pdf', width=10, height=20)
ggplot(melt(test.transformed, id.vars = "Shifting_Prom"),
       aes(x = variable, y = value, col = Shifting_Prom)) +
  geom_violin() + facet_wrap( ~variable, ncol = 2, scales = "free") +
  theme_bw()
dev.off()
######################################################
## compact
library(dplyr)
all.peaks.pgc1 <- all.peaks.pgc %>% group_by(geneId) %>%
  summarise(count = n()) %>% data.frame()

## get coordinates of peaks (esclude promoters)
## assing gene by distance
## merge with rnaseq


### table

# ENSDARG00000104632  GeneType  TPM #enhancers  

###################### plot

# plot n TATA shift correlation
library(dplyr)
library(reshape2)
library(ggplot2)

## no ATAC
ggplot(final_data_3_list_no_Atac[[2]], aes(x = log(Expr_var_PGC), y = width, fill = Shift_Prim5,
                            colour = Shift_Prim5)) +
  geom_point() + theme_bw()

## ATAC
ggplot(final_data_3_list_Atac[[3]], aes(x = Expr_var_PGC, y = num_ATAC_peaks, fill = Shift_Prim5,
                                           colour = Shift_Prim5)) +
  geom_point() + theme_bw()

# 

for(x in range(1:length(final_data_3_list))){
  #final_data_3_list[[x]][is.na(final_data_3_list[[x]])] <- 0
  m <- final_data_3_list[[x]][,c(17:22)]
  mm <- melt(m)
  f <- data.frame(Shifting_Scores = final_data_3_list[[x]][,33])
  # append shifting scores
  a <- b <- c <- d <- e <- f
  sh_sc <- bind_rows(a,b,c,d,e,f)
  
  mm$Shifting_Scores <- sh_sc$Shifting_Scores
  
  #  dot plot
  my_name <- names(final_data_3_list)[x]
  #pdf(paste0('~/Desktop/Postdoc/Data_Analysis/ML_PGC_promoters/Plots/Corr_TATApos_ShiftScore', my_name, '.png'))
  ggplot(mm, aes(x = value, y = Shifting_Scores, fill = variable, color = variable)) + 
    geom_point() + xlim(1,150) + ylim(0.59,1) + theme_bw()
  #dev.off()
}

## exctract shifting prom with tata at -120

#####
library(viridis)
vallie <- sample(viridis(32))

# PLOT: number of ATAC peaks per gene
final_data_plot <- final_data_TBP_5 %>%
group_by(GeneType) %>%
add_tally(name = "sample_n")  %>%
group_by(Num_ATAC_peaks, GeneType) %>%
summarise(N = n(), prop = (N/sample_n[1])*100,
n_sample = sample_n[1])

ggplot(final_data_plot, aes(x = GeneType, y = prop, fill = as.character(Num_ATAC_peaks),
                            colour = as.character(Num_ATAC_peaks))) +
  geom_bar(stat="identity") + scale_fill_manual(values = vallie) +
  scale_color_manual(values = vallie) + theme_bw()

# PLOT: number of 
final_data_plot <- final_data_TBP_5 %>%
group_by(GeneType) %>%
summarise(Mean_distance = mean(Distance_ATAC_Mean))

ggplot(final_data_plot, aes(x = GeneType, y = Mean_distance)) +
geom_bar(stat="identity") + scale_fill_manual(values = vallie) +
scale_color_manual(values = vallie) + theme_bw()

# PLOT: proportion of gene biotype
final_data_plot <- final_data_TBP_5 %>%
group_by(GeneType) %>%
add_tally(name = "sample_n")  %>%
group_by(gene_biotype, GeneType) %>%
summarise(N = n(), prop = (N/sample_n[1])*100,
n_sample = sample_n[1])

final_data_plot <- final_data_plot[-c(13,14,15),]
ggplot(final_data_plot, aes(x = GeneType, y = prop, fill = as.character(gene_biotype),
                            colour = as.character(gene_biotype))) +
  geom_bar(stat="identity") + scale_fill_manual(values = vallie) +
  scale_color_manual(values = vallie) + theme_bw()

final_data_plot[is.na(final_data_plot)] <- "No promoter"

# PLOT: num shifting promoters
final_data_plot <- final_data_TBP_5 %>%
group_by(GeneType) %>%
add_tally(name = "sample_n")  %>%
group_by(Shifting_Prom, GeneType) %>%
summarise(N = n(), prop = (N/sample_n[1])*100,
n_sample = sample_n[1])
final_data_plot[is.na(final_data_plot)] <- "No promoter"

ggplot(final_data_plot, aes(x = GeneType, y = prop, fill = as.character(Shifting_Prom),
                            colour = as.character(Shifting_Prom))) +
  geom_bar(stat="identity") + scale_fill_manual(values = vallie) +
  scale_color_manual(values = vallie) + theme_bw()

# PLOT num broad/sharp promoters
final_data_plot <- final_data_TBP_5 %>%
  group_by(GeneType) %>%
  add_tally(name = "sample_n")  %>%
  group_by(Width_Class, GeneType) %>%
  summarise(N = n(), prop = (N/sample_n[1])*100,
            n_sample = sample_n[1])
final_data_plot[is.na(final_data_plot)] <- "No promoter"

ggplot(final_data_plot, aes(x = GeneType, y = prop, fill = as.character(Width_Class),
                            colour = as.character(Width_Class))) +
  geom_bar(stat="identity") + scale_fill_manual(values = vallie) +
  scale_color_manual(values = vallie) + theme_bw()

# PLOT
final_data_plot <- final_data_TBP_5 %>%
  group_by(GeneType)

ggplot(final_data_TBP_5, aes(x = interquantile_width, fill = GeneType)) + geom_density(alpha = 0.4) +
  theme_bw(base_size=2)

