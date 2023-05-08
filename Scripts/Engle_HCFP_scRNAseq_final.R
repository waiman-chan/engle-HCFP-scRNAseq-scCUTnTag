
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)



##Basic info##
run_name <- "cfw_cre1dup_e95-125_DUP_FilterthenMerge_Hoxb1-IsL_PC40_C7mm39"

data_name <- c("cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td", "cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td",
               "cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td", "cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td",
               "cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td", "cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td",
               "cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td", "cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td")

sample_name <- c("e095_wt", "e095_dup",
                 "e105_wt", "e105_dup",
                 "e115_wt", "e115_dup",
                 "e125_wt", "e125_dup")

sample_name <- c("e095_wt", "e105_wt",
                 "e115_wt", "e125_wt")

sample_name <- c("e095_dup", "e105_dup",
                 "e115_dup", "e125_dup")

####1: Load data and create Seurat object ######################################
for (file in data_name){ 
  seurat_data <- Read10X(data.dir = paste0("data/cfw/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

##Add sample column
cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td[["sample"]] <- 'e095_wt'
cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td[["sample"]] <- 'e105_wt'
cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td[["sample"]] <- 'e115_wt'
cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td[["sample"]] <- 'e125_wt'
cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td[["sample"]] <- 'e095_dup'
cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td[["sample"]] <- 'e105_dup'
cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td[["sample"]] <- 'e115_dup'
cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td[["sample"]] <- 'e125_dup'


####2: FILTER CELLS USING QUALITY MATRICS##############################
## FILTER CELLS USING QUALITY MATRICS

# Compute percent mito ratio
cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td$mitoRatio <- PercentageFeatureSet(object = cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td, pattern = "^mt-") 
cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td$mitoRatio <- cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td@meta.data$mitoRatio / 100
  
cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td$mitoRatio <- PercentageFeatureSet(object = cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td, pattern = "^mt-") 
cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td$mitoRatio <- cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td@meta.data$mitoRatio / 100

cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td$mitoRatio <- PercentageFeatureSet(object = cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td, pattern = "^mt-") 
cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td$mitoRatio <- cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td@meta.data$mitoRatio / 100

cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td$mitoRatio <- PercentageFeatureSet(object = cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td, pattern = "^mt-") 
cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td$mitoRatio <- cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td@meta.data$mitoRatio / 100

cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td$mitoRatio <- PercentageFeatureSet(object = cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td, pattern = "^mt-")
cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td$mitoRatio <- cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td@meta.data$mitoRatio / 100

cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td$mitoRatio <- PercentageFeatureSet(object = cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td, pattern = "^mt-")
cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td$mitoRatio <- cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td@meta.data$mitoRatio / 100

cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td$mitoRatio <- PercentageFeatureSet(object = cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td, pattern = "^mt-")
cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td$mitoRatio <- cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td@meta.data$mitoRatio / 100

cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td$mitoRatio <- PercentageFeatureSet(object = cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td, pattern = "^mt-")
cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td$mitoRatio <- cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td@meta.data$mitoRatio / 100


#filter out mito
cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td <- subset(x = cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td, subset = mitoRatio < 0.05)
cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td <- subset(x = cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td, subset = mitoRatio < 0.05)
cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td <- subset(x = cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td, subset = mitoRatio < 0.05)
cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td<- subset(x = cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td, subset = mitoRatio < 0.05)
cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td <- subset(x = cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td, subset = mitoRatio < 0.05)
cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td <- subset(x = cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td, subset = mitoRatio < 0.05)
cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td <- subset(x = cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td, subset = mitoRatio < 0.05)
cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td<- subset(x = cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td, subset = mitoRatio < 0.05)


#count distribution - e9.5 wt
jpeg(paste0('E09.5_wt','violinplot_count-feature-mito.jpg'))
VlnPlot(cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td, features = c('nCount_RNA', 'nFeature_RNA', 'mitoRatio'), pt.size = 0)
dev.off()

jpeg(paste0('E09.5_wt','scatter_count-feature.jpg'))
FeatureScatter(
  cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

q_umi_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td$nCount_RNA, probs = c(.01, .05, .95, .99))
cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td <- subset(cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td, subset = nCount_RNA <= q_umi_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td[4] & 
                                          nCount_RNA >= q_umi_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td[1]) #top and bottom 1%

q_feature_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td$nFeature_RNA, probs = c(.01, .05, .95, .99))
q_feature_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td

#count distribution - e9.5 dup
jpeg(paste0('E09.5_dup','violinplot_count-feature-mito.jpg'))
VlnPlot(cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td, features = c('nCount_RNA', 'nFeature_RNA', 'mitoRatio'), pt.size = 0)
dev.off()

jpeg(paste0('E09.5_dup','scatter_count-feature.jpg'))
FeatureScatter(
  cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

q_umi_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td$nCount_RNA, probs = c(.01, .05, .95, .99))
cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td <- subset(cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td, subset = nCount_RNA <= q_umi_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td[4] &
                                          nCount_RNA >= q_umi_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td[1]) #top and bottom 1%

q_feature_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td$nFeature_RNA, probs = c(.01, .05, .95, .99))
q_feature_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td


#count distribution - e10.5 wt
jpeg(paste0('E10.5_wt','violinplot_count-feature-mito.jpg'))
VlnPlot(cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td , features = c('nCount_RNA', 'nFeature_RNA', 'mitoRatio'), pt.size = 0)
dev.off()

jpeg(paste0('E10.5_wt','scatter_count-feature.jpg'))
FeatureScatter(
  cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td , "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

q_umi_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td <-quantile(cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td$nCount_RNA, probs = c(.01, .05, .95, .99))
cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td  <- subset(cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td, subset = nCount_RNA <= q_umi_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td[4] & 
                                            nCount_RNA >= q_umi_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td[1]) #top and bottom 1%

q_feature_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td$nFeature_RNA, probs = c(.01, .05, .95, .99))
q_feature_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td


#count distribution - e10.5 dup
jpeg(paste0('E10.5_dup','violinplot_count-feature-mito.jpg'))
VlnPlot(cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td , features = c('nCount_RNA', 'nFeature_RNA', 'mitoRatio'), pt.size = 0)
dev.off()

jpeg(paste0('E10.5_dup','scatter_count-feature.jpg'))
FeatureScatter(
  cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td , "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

q_umi_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td <-quantile(cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td$nCount_RNA, probs = c(.01, .05, .95, .99))
cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td  <- subset(cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td, subset = nCount_RNA <= q_umi_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td[4] &
                                            nCount_RNA >= q_umi_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td[1]) #top and bottom 1%

q_feature_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td$nFeature_RNA, probs = c(.01, .05, .95, .99))
q_feature_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td


#count distribution - e11.5 wt
jpeg(paste0('E11.5_wt','violinplot_count-feature-mito.jpg'))
VlnPlot(cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td, features = c('nCount_RNA', 'nFeature_RNA', 'mitoRatio'), pt.size = 0)
dev.off()

jpeg(paste0('E11.5_wt','scatter_count-feature.jpg'))
FeatureScatter(
  cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

q_umi_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td$nCount_RNA, probs = c(.01, .05, .95, .99))
cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td <- subset(cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td, subset = nCount_RNA < q_umi_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td[4] & 
                                                  nCount_RNA > q_umi_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td[1]) #top and bottom 1%

q_feature_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td$nFeature_RNA, probs = c(.01, .05, .95, .99))
q_feature_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td 

#count distribution - e11.5 dup
jpeg(paste0('E11.5_dup','violinplot_count-feature-mito.jpg'))
VlnPlot(cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td, features = c('nCount_RNA', 'nFeature_RNA', 'mitoRatio'), pt.size = 0)
dev.off()

jpeg(paste0('E11.5_dup','scatter_count-feature.jpg'))
FeatureScatter(
  cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

q_umi_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td$nCount_RNA, probs = c(.01, .05, .95, .99))
cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td <- subset(cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td, subset = nCount_RNA < q_umi_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td[4] &
                                                  nCount_RNA > q_umi_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td[1]) #top and bottom 1%

q_feature_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td$nFeature_RNA, probs = c(.01, .05, .95, .99))
q_feature_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td

#count distribution - e12.5 wt
jpeg(paste0('E12.5_wt','violinplot_count-feature-mito.jpg'))
VlnPlot(cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td, features = c('nCount_RNA', 'nFeature_RNA', 'mitoRatio'), pt.size = 0)
dev.off()

jpeg(paste0('E12.5_wt','scatter_count-feature.jpg'))
FeatureScatter(
  cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

q_umi_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td$nCount_RNA, probs = c(.01, .05, .95, .99))
cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td <- subset(cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td, subset = nCount_RNA <= q_umi_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td[4] & 
                                                  nCount_RNA >= q_umi_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td[1]) #top and bottom 1%

q_feature_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td$nFeature_RNA, probs = c(.01, .05, .95, .99))
q_feature_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td 

#count distribution - e12.5 dup
jpeg(paste0('E12.5_dup','violinplot_count-feature-mito.jpg'))
VlnPlot(cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td, features = c('nCount_RNA', 'nFeature_RNA', 'mitoRatio'), pt.size = 0)
dev.off()

jpeg(paste0('E12.5_dup','scatter_count-feature.jpg'))
FeatureScatter(
  cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

q_umi_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td$nCount_RNA, probs = c(.01, .05, .95, .99))
cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td <- subset(cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td, subset = nCount_RNA <= q_umi_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td[4] &
                                                  nCount_RNA >= q_umi_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td[1]) #top and bottom 1%

q_feature_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td <- quantile(cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td$nFeature_RNA, probs = c(.01, .05, .95, .99))
q_feature_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td

####2. Merge objects############################################################
merged_seurat <- merge(cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td,
                       y = c(cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td,
                             cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td, cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td,
                             cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td, cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td,
                             cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td, cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td),
                       add.cell.id = sample_name)
# 
merged_seurat <- merge(cfw_cre1dup_e95_wt_2_mm39c7_excluded_gfp_td,
                       y = c(cfw_cre1dup_e105_wt_3_mm39c7_excluded_gfp_td, cfw_cre1dup_e115_wt_combined_mm39c7_excluded_gfp_td,
                             cfw_cre1dup_e125_wt_combined_mm39c7_excluded_gfp_td),
                       add.cell.id = sample_name)
# # # 
merged_seurat <- merge(cfw_cre1dup_e95_dup_2_mm39c7_excluded_gfp_td,
                       y = c(cfw_cre1dup_e105_dup_3_mm39c7_excluded_gfp_td, cfw_cre1dup_e115_dup_combined_mm39c7_excluded_gfp_td,
                             cfw_cre1dup_e125_dup_combined_mm39c7_excluded_gfp_td),
                       add.cell.id = sample_name)

####7. OPTIONAL STEPS: Subsetting by gene expression - added by Jo#############
#7. OPTIONAL STEPS: Subsetting by gene expression - added by Jo

#### filtering cells by gene expression
filtered_seurat <- merged_seurat

filtered_seurat_subset <- subset(x = filtered_seurat, subset = Hoxb1 > 0 | Isl1 >0)

for (s in sample_name){
  n_cells <- subset(x = filtered_seurat_subset, subset = (sample == s))
  print(s)
  print(ncol(assign(paste0(s,"_n_cells_Hoxb1-Isl"), n_cells)))
}


####8. Explore sources of unwanted variation##################################
#8. Explore sources of unwanted variation
# do this if filter then merge and for Hoxb1 and IsL cells only
filtered_seurat <- filtered_seurat_subset 
seurat_phase <- NormalizeData(filtered_seurat)

##8a. Evaluating effects of cell cycle
# Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle. Read it in with:

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 
cell_cycle_genes <- read.csv(text = cc_file)

#------------------------------------------------------------------#
# CONVERT ENSEMBL IDS TO GENE NAME
# Install AnnotationHub  
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("AnnotationHub")

library(AnnotationHub)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, 
                                       by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
#View(seurat_phase@meta.data)   

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, #replaced by SCTransform
                                     selection.method = "vst",
                                     nfeatures = 2000, # Tommy uses 2500
                                     verbose = FALSE)

## check number of variable genes
length(seurat_phase@assays$RNA@var.features)

# Identify the 10 most highly variable genes
top10 <- head(seurat_phase@assays$RNA@var.features, 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_phase)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
jpeg(paste0(run_name, '_PCA_byPhase.jpg'))
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()

##8b. Evaluating effects of mitochondrial expression

# Check quantile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))

#plot the PCA similar to how we did with cell cycle regression
# Plot the PCA colored by mitoFr
jpeg(paste0(run_name, '_byMito.jpg'))
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")
dev.off()
####9. Normalization and regressing out sources of unwanted variation using SCTransform####
##Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()

#9. Normalization and regressing out sources of unwanted variation using SCTransform
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat_no_regression <- SplitObject(seurat_phase, split.by = "sample")

split_seurat_no_regression <- split_seurat_no_regression[sample_name]

# Adjusting memory
#options(future.globals.maxSize = 4000 * 1024^2)

#perform the sctransform on all samples
for (i in 1:length(split_seurat_no_regression)) {
  split_seurat_no_regression[[i]] <- SCTransform(split_seurat_no_regression[[i]])
}

# Save the split seurat object
saveRDS(split_seurat_no_regression, paste0(run_name, "_split_seurat_no_regression.rds"))

####11. Integration##########################################################
#11. Integration

##11a: Select the most variable features to use for integration
integ_features_no_regression <- SelectIntegrationFeatures(object.list = split_seurat_no_regression, 
                                                          nfeatures = 3000)

##11b: Prepare the SCT list object for integration
split_seurat_no_regression <- PrepSCTIntegration(object.list = split_seurat_no_regression, 
                                                 anchor.features = integ_features_no_regression)

##11c: Find best buddies - can take a while to run
integ_anchors_no_regression <- FindIntegrationAnchors(object.list = split_seurat_no_regression, 
                                                      normalization.method = "SCT", 
                                                      anchor.features = integ_features_no_regression)

##11d: Integrate across conditions
seurat_integrated_no_regression <- IntegrateData(anchorset = integ_anchors_no_regression, 
                                                 normalization.method = "SCT")

##12e: Save integrated seurat object
saveRDS(seurat_integrated_no_regression, 
        paste0(run_name, "_integrated_seurat_no_regression.rds"))

####12. UMAP visualization###################################################
#12. UMAP visualization

##12a: Run PCA
seurat_integrated_no_regression <- RunPCA(object = seurat_integrated_no_regression)

##12b: Plot PCA
jpeg(paste0(run_name, '_PCA_separate.jpg'))
PCAPlot(seurat_integrated_no_regression,
        split.by = "sample")  
dev.off()

##12c: Run UMAP
seurat_integrated_no_regression <- RunUMAP(seurat_integrated_no_regression, 
                                           dims = 1:40,
                                           reduction = "pca")

##12d: Plot UMAP   
jpeg(paste0(run_name, '_UMAP_overlay.jpg'))
DimPlot(seurat_integrated_no_regression) 
dev.off()

# Plot UMAP split by sample
jpeg(paste0(run_name, '_UMAP_separate.jpg'))
DimPlot(seurat_integrated_no_regression,
        split.by = "sample",
        ncol = 2) 
dev.off()

##12e: Save integrated seurat object
saveRDS(seurat_integrated_no_regression, 
        paste0(run_name, "_integrated_seurat_no_regression.rds"))




#### SKIPPED Step13 ####
# Plot the elbow plot
jpeg(paste0(run_name,'_elbow.jpg'))
ElbowPlot(object = seurat_integrated_no_regression, 
          ndims = 40)
dev.off()
####14. Clustering cells###################################################
#14. Clustering cells

##14a: Determine the K-nearest neighbor graph
seurat_integrated_no_regression <- FindNeighbors(object = seurat_integrated_no_regression, 
                                                 dims = 1:40)

##14b: Determine the clusters for various resolutions                                
seurat_integrated_no_regression <- FindClusters(object = seurat_integrated_no_regression,
                                                resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.4))

##14c: Explore resolutions
seurat_integrated_no_regression@meta.data %>% 
  View()

for (i in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.4)) {
  jpeg(paste0(run_name,"_UMAP_", i, ".jpg"))
  seurat_integrated <- FindClusters(seurat_integrated_no_regression, resolution = i)
  print(DimPlot(seurat_integrated,reduction = 'umap') + ggtitle('Resolution:', i))
  dev.off()
}

##14d: Assign identity of clusters
Idents(object = seurat_integrated_no_regression) <- "integrated_snn_res.0.2" 
# need to change to desired resolution

##14e: Plot the final UMAP
# UMAP of cells in each cluster by sample
jpeg(paste0(run_name,'_UMAP_bySample.jpg'))
DimPlot(seurat_integrated_no_regression, 
        label = TRUE, 
        split.by = "sample",
        ncol = 2)  + NoLegend()
dev.off()

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated_no_regression, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)
write.csv(n_cells, file = paste0(run_name, '_ncells.csv'))

#-------------------------------------------------------------#
# Segregation of clusters by various sources of uninteresting variation

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

jpeg(paste0(run_name, '_feature_byMetrics.jpg'))
FeaturePlot(seurat_integrated_no_regression, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()

#-------------------------------------------------------------#
# Exploration of the PCs driving the different clusters
# Not included here

# Save integrated seurat object
saveRDS(seurat_integrated_no_regression, 
        paste0(run_name,"_integrated_seurat_no_regression.rds"))
