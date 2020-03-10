\name{Output}
\alias{Output}
\title{Take the bin number of the detected TAD boundaries as input, and output the
       TAD boundaries in an optional form between text file and graphical heatmap}
\usage{
Output(df_result, species, chr, resolution, outxtfile ="./result", bheatmap=F, heatmapfile=NULL, hicmat=NULL)

df_result       a dataframe consisting of the bin number of the detected TAD
                boundaries by TADBD function.
species         species name, including "hg19", "hg38", "mm9" and "mm10".
chr             chromosome name.
resolution      resolution of of Hi-C matrix.
outxtfile       the path and name of output text file.
bheatmap        if a heatmap with TAD tracks is shown, bheatmap = T, otherwise
                bheatmap = F(default).
heatmapfile     the path and name of output heatmap file, the argument needs to be
                specifie only when bheatmap = T.
hicmat          a dense contact matrix by DataLoad, the argument needs to be
                specifie only when bheatmap = T.
}
\description{Take the bin number of the detected TAD boundaries as input, and output the TAD boundaries in an optional form between text file and graphical heatmap.
}
\examples{
#Clear all objects from the current workspace
rm(list=ls())
#Load R package TADBD
library(TADBD)
#Configuration of the parameters, including species, chromsome and resolution
species <- "hg19"
chr <- "chr18"
resolution <- 50000
#Close scientific notation
options(scipen = 999)
#Specify Hi-C data to be loaded
data(hicdata)
#Load a Hi-C contact matrix file in a dense format
hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)
#Detect TAD boundaries on the loaded contact matrix using TADBD method
df_result <- TADBD(hicmat)
#Output two text files, one is for detected TAD boundaries, the other for intermediate peaks 
Output(df_result, species, chr, resolution, outxtfile="./result")
#Output a text file with TAD boundary coordinates, as well as a heatmap with TAD boundary tracks
Output(df_result, species, chr, resolution, outxtfile="./result", bheatmap = F, heatmapfile="./heatmap", hicmat)
}
