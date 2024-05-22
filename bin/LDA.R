library("tidyverse")
library("ggplot2")
library("reshape2")
library("dplyr")
library("MASS")

##################################################################################################
# Workflow for data transformation and Linear Discriminant Analysis

output_path <- '../Output'
final_data_5_list <- readRDS('data/ML_Final_Complete_List_PGC_Soma_Wt.rds')

## prepare the data for ML
.data_clean <- function(x){
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

r <- lapply(final_data_5_list, .data_clean)

## plot correlations
r1 <- r[[2]]
r1$Width_Class <- NULL

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

# PLOT: Mean distance
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

