# OVERVIEW 
The repository includes the final scripts and sample data set used for analyses that are associated with the manuscript "Non-coding variants alter GATA2 expression in rhombomere 4 motor neurons and cause dominant hereditary congenital facial paresis" by Tenney et al. 

Experiments include
1. scRNAseq data analysis for e9.5, e10.5, e11.5, and e12.5 wild-type and cRE1dup/+ mutant samples.
2. scCUT&Tag data analysis for e10.5 wild-type and Fam5 snv/snv mutant samples

# SAMPLE DATA SET

The folder "Sample data" contains:
1. A set of wild-type and cRE1dup/+ mutant E12.5 data files for scRNAseq data analysis
2. A set of wild-type and Fam5snv/snv mutant E10.5 data files for scCUT&Tag analysis

# COMPLETE DATA SET
A complete set of our Single-cell RNA and CUT&Tag sequencing data are available through NCBI GEO SuperSeries record GSE223274. 

#CUSTOM CODE
The folder "Scripts" contains:
1. The final code used in our scRNAseq data analysis
2. The final code used in our scCUT&Tag data analysis

# HARDWARE REQUIREMENT

Scripts provided here require only a standard computer with enough RAM to support all the operations. 
Non-standard hardware is not required.

# SOFTWARE REQUIREMENT

Scripts can be run on R alone or within R Studio. The scripts have been tested on the following version of software: 
R (version 4.2.1, v4.2.11.2)
R Studio (version 2022.12.0+353)

# INSTALLATION GUIDE
## Software

A copy of R, and R Studio can be obtained and installed here: https://posit.co/download/rstudio-desktop/. Typical installation time of both R and R studio on a normal computer should not exceed 30 minutes. 

No other software is required for running the scripts. 

## R libraries
Please be sure to download and load all the necessary packages before running the script. The libraries needed are listed at the top of each script file. 

# INSTRUCTIONS/DESCRIPTION
## scRNAseq
1. Create a project name as a variable. This will be added to all the file names generated in this analysis.
2. Create a vector of your data file name. This is to create Seurat objects.
3. Create a vector of sample names corresponding to the data file name. This is to label the Seurat objects.
4. Load data and create Seurat Object: enter the file path to your data files and use Read10x and CreateSeuratObject to create Seurat object
5. For each object, add a column "sample" in the metadata of the object and add sample name to each row
6. For each object, filter cell suing quality metrics: percent of mitrochrondia reads; review count distribution and remove top and bottom 1% of cells. 
7. Merge objects.
8. Subset cells with Hoxb1 and IsL expression and save as a filtered object. 
9. Normalize data
10. Score cells for cell cycle.
11. Find Variable Features
12. Scale data and run PCA.
13. Evaluate effects of mitochondrial expression
14. Plot PCA by mito expression.
15. SCTransform.
16. Save object (.rds file)
17. Data integration.
18.Save object as an integrated object (.rds file)
19. UMap utilization
20. Save object as an integrated object (.rds file). This is the final object used for analysis in BBrowser3 (BioTuring).
To execute all code in the script using our full data set, it should not take more than an hour. 

## scCUT&Tag
This script will allow to generate the pseudo-bulk figure showed in the paper.
1. Import the peak files for all samples.
2. Convert the peaks to genomic ranges.
3. Generate a combined peak files and filter based to the peak width [>20 and <10,000]
4. Generate the fragment object by importing the metadata file (singlecell.csv file) and filter low count cells.
5. Quantify peaks in each individual dataset using FeatureMatrix function.
6. Create a Seurat object with all the previously generated files.
7. Create a merged Seurat object with all replicates (Remember to add a column with the dataset definition before merging).
8. Run dimensionality reduction and TFIDF function using the Signac function RunTFIDF, FindTopFeature, RunSVD and RunUMAP (with dims=2:50 and reduction=”lsi”).
9. Plot the combined dataset in the UMAP space.
10. Add genomic coordinate and annotation to the dataset.
11. Generate a pseudo bulk Coverage plot of the region of the interest (mm10, chr6-88186706-88250381).
12. Save and export the Coverage plot figure.

# Zenodo Link
https://zenodo.org/badge/latestdoi/635023782
