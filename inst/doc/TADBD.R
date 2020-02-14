## ----set-options, echo=FALSE, cache=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 400)

## ---- eval = FALSE, message=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  # if (!requireNamespace("BiocManager", quietly=TRUE))
#  #     install.packages("BiocManager")
#  # BiocManager::install("TADBD")
#  #devtools::install_github("bioinfo-lab/TADBD")
#  #library(TADBD)

## ---- echo = FALSE, warning = FALSE, message = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#rm(list=ls())
#library(TADBD)
#species <- "hg19"
#chr <- "chr18"
#resolution <- 50000
#options(scipen = 999)
#data(hicdata)
#hicmat <- DataLoad(hicdata, bsparse = T, species, chr, resolution)

## ---- echo = FALSE, warning = FALSE, message = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(TADBD)
species <- "hg19"
chr <- "chr18"
resolution <- 50000
options(scipen = 999)
data(hicdata)
hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)

## ---- message = FALSE, warning = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(TADBD)
#specify the species, resolution, and chromosome
species <- "hg19"
chr <- "chr18"
resolution <- 50000
options(scipen = 999)
data(hicdata)
#Transfer the input format into an acceptable format
hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)
#Output the bin number of TAD boundaries on the contact matrix 
df_result <- TADBD(hicmat)

## ---- message = FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(TADBD)
species <- "hg19"
chr <- "chr18"
resolution <- 50000
options(scipen = 999)
data(hicdata)
hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)
df_result <- TADBD(hicmat)
#Output coordinates of TAD boundary 
Output(df_result, species, chr, resolution, outxtfile="./result")


## ---- message = FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(TADBD)
species <- "hg19"
chr <- "chr18"
resolution <- 50000
options(scipen = 999)
data(hicdata)
hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)
df_result <- TADBD(hicmat)
#Output heatmap and TAD boundary tracks
Output(df_result, species, chr, resolution, outxtfile="./result", bheatmap = T, heatmapfile="./heatmap", hicmat)

