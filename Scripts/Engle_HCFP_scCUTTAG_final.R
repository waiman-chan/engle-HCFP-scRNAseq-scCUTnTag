#01/20/2023
#Rversion 4.1.2
#This scripts have been used to generate the scCUT&TAG plot for the analysis
#Author:Alessandro Di Gioia (alessandro.digioia2@gmail.com). For any additional information contact the author

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
set.seed(1234)


#Read peak set for both samples

peaks.WT_1 <- read.table("../../../sdigioia/Desktop/BCH_work/GFP-NR2F1_2/peaks.bed", col.names = c("chr", "start", "end"))
peaks.SNV_1 <- read.table("../../../sdigioia/Desktop/BCH_work/CGQ_SNV_Replicate1/peaks.bed",col.names = c("chr", "start", "end"))
peaks.WT_2 <- read.table("../../../sdigioia/Desktop/BCH_work/GFP_NR2F1_1/peaks.bed", col.names = c("chr", "start", "end"))
peaks.SNV_2 <- read.table("../../../sdigioia/Desktop/BCH_work/CGQ_SNV_Replicate2/peaks.bed",col.names = c("chr", "start", "end"))

#convert to genomic range

gr.WT_1 <- makeGRangesFromDataFrame(peaks.WT_1)
gr.SNV_1 <- makeGRangesFromDataFrame(peaks.SNV_1)
gr.WT_2 <- makeGRangesFromDataFrame(peaks.WT_2)
gr.SNV_2 <- makeGRangesFromDataFrame(peaks.SNV_2)

#Create unified set of peaks to quantify in each dataset

combined.peaks <- reduce(x=c(gr.WT_1,gr.SNV_1, gr.WT_2, gr.SNV_2))
#Filter bad peaks based on lenght

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
combined.peaks

#Create the fragment object

#Load metadata

md.WT_1 <- read.table(
  file = "../../../sdigioia/Desktop/BCH_work/GFP-NR2F1_2/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1)[-1,]


md.SNV_1 <- read.table(
  file = "../../../sdigioia/Desktop/BCH_work/CGQ_SNV_Replicate1/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1)[-1,] 

md.WT_2 <- read.table(
  file = "../../../sdigioia/Desktop/BCH_work/GFP_NR2F1_1/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1)[-1,]

md.SNV_2 <- read.table(
  file = "../../../sdigioia/Desktop/BCH_work/CGQ_SNV_Replicate2/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1)[-1,] 

#Perform filtering on low count cells
md.WT_1 <- md.WT_1[md.WT_1$passed_filters > 1000,]
md.SNV_1 <- md.SNV_1[md.SNV_1$passed_filters > 1000,]
md.WT_2 <- md.WT_2[md.WT_2$passed_filters > 1000,]
md.SNV_2 <- md.SNV_2[md.SNV_2$passed_filters > 1000,]


#Create the fragment objects
#Note that I got an error for the cell count let's try using the 
frags.WT_1 <- CreateFragmentObject(
  path = "../../../sdigioia/Desktop/BCH_work/GFP-NR2F1_2/fragments.tsv.gz",
)


frags.SNV_1 <- CreateFragmentObject(
    path = "../../../sdigioia/Desktop/BCH_work/CGQ_SNV_Replicate1/fragments.tsv.gz",
  )

frags.WT_2 <- CreateFragmentObject(
  path = "../../../sdigioia/Desktop/BCH_work/GFP_NR2F1_1/fragments.tsv.gz",
)


frags.SNV_2 <- CreateFragmentObject(
  path = "../../../sdigioia/Desktop/BCH_work/CGQ_SNV_Replicate2/fragments.tsv.gz",
)

#Quantify peaks in each dataset
WT.counts_1 <- FeatureMatrix(
  fragments = frags.WT_1, 
  features = combined.peaks
)


SNV.counts_1 <- FeatureMatrix(
  fragments = frags.SNV_1, 
  features = combined.peaks
)

WT.counts_2 <- FeatureMatrix(
  fragments = frags.WT_2, 
  features = combined.peaks
)


SNV.counts_2 <- FeatureMatrix(
  fragments = frags.SNV_2, 
  features = combined.peaks
)

#Create the object

WT_assay_1 <- CreateChromatinAssay(
  WT.counts_1, 
  fragments = frags.WT_1)

WT_1 <- CreateSeuratObject(WT_assay_1, assay ="scCutTAG", meta.data = md.WT_1)



SNV_assay_1 <- CreateChromatinAssay(SNV.counts_1, fragments = frags.SNV_1)
SNV_1 <- CreateSeuratObject(SNV_assay_1, assay = "scCutTAG", meta.data = md.SNV_1)

WT_assay_2 <- CreateChromatinAssay(
  WT.counts_2, 
  fragments = frags.WT_2)

WT_2 <- CreateSeuratObject(WT_assay_2, assay ="scCutTAG", meta.data = md.WT_2)


SNV_assay_2 <- CreateChromatinAssay(SNV.counts_2, fragments = frags.SNV_2)
SNV_2 <- CreateSeuratObject(SNV_assay_2, assay = "scCutTAG", meta.data = md.SNV_2)

######Merge the objects#####

# add information to identify the dataset

WT_1$dataset <- 'WT'
SNV_1$dataset <- 'SNV'
WT_2$dataset <- 'WT'
SNV_2$dataset <- 'SNV'


# add information on replicate
WT_1$replicate <- 'WT_rep1'
SNV_1$replicate <- 'snv_rep1'
WT_2$replicate <- 'WT_rep2'
SNV_2$replicate <- 'snv_rep2'


#merge all dataset adding cell ID to make sure that cell names are unique

combined <-  merge(
  x = WT_1,
  y=list(WT_2,SNV_1,SNV_2),
  add.cell.id = c("WT_cell_rep1","WT_cell_rep2", "SNV_cell_rep1", "SNV_cell_rep2")
)
combined[["scCutTAG"]]


combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
DimPlot(combined, group.by = 'replicate')

#add annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to mm10
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(combined) <- annotations

##Check the specific region
CoveragePlot(
  object = combined,
  region = c("chr6-88186706-88250381"),
  group.by= 'dataset'
)

#plot by replicate
CoveragePlot(
  object = combined,
  region = c("chr6-88186706-88250381"),
  group.by= 'replicate'
)
##########################################################################################


