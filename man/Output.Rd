\name{Output}
\alias{Output}
\title{Take the bin number of the detected TAD boundaries as input, and output the
       TAD boundaries in an optional form between text file and graphical heatmap}
\usage{
Output(df_result, species, chr, resolution, outxtfile ="./result", bheatmap=FALSE, heatmapfile=NULL, hicmat=NULL, map_start = 0, map_end = 10000000, l_color="blue", l_width=2.5)

df_result       a dataframe consisting of the bin number of the detected TAD
                boundaries by TADBD function.
species         species name, including "hg19", "hg38", "mm9" and "mm10".
chr             chromosome name.
resolution      resolution of of Hi-C matrix.
outxtfile       the path and name of output text file.
bheatmap        if a heatmap with TAD tracks is shown, bheatmap = TRUE, otherwise
                bheatmap = FALSE(default).
heatmapfile     the path and name of output heatmap file, the argument needs to be
                specifie only when bheatmap = TRUE.
hicmat          a dense contact matrix by DataLoad, the argument needs to be
                specifie only when bheatmap = TRUE.
map_start       the starting coordinate of the heatmap,value is an integer multiple of resolution and the default is 0.
map_end         the ending coordinate of the heatmap,value is an integer multiple of resolution and the default is 10000000.
l_color         color of the TAD track.
l_width         width of the TAD track.
}
\description{Take the bin number of the detected TAD boundaries as input, and output the TAD boundaries in an optional form between text file and graphical heatmap.
}
\examples{
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
hicmat <- DataLoad(hicdata, bsparse = FALSE, species, chr, resolution)
#Detect TAD boundaries on the loaded contact matrix using TADBD method using default parameter configuration, that is template.sizes = c(4,5,6), bstatfilter = TRUE
df_result <- TADBD(hicmat)
#Output two text files, one is for detected TAD boundaries, the other for intermediate peaks
Output(df_result, species, chr, resolution, outxtfile="./result")
#Output two text files and a heatmap with TAD boundary tracks, the parameters of heatmap include starting and ending coordinates, as well as the color and the width of tracks
Output(df_result, species, chr, resolution, outxtfile="./result", bheatmap = TRUE, heatmapfile="./heatmap", hicmat, map_start=0, map_end=10000000, l_color="bule", l_width=2.5)
}
