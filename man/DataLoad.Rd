\name{DataLoad}
\alias{DataLoad}
\title{Load a Hi-C matrix file in sparse or dense format, and output the dense
       contact matrix}
\usage{
DataLoad(hicdata, bsparse, species, chr, resolution)}
\arguments{
\item{hicdata}{a Hi-C matrix file.}
\item{bsparse}{if the Hi-C matrix is in sparse format, bsparse = TRUE, otherwise bsparse = FALSE.}
\item{species}{species name, including "hg19", "hg38", "mm9" and "mm10", the argument needs to be specified only when bsparse = TRUE.}
\item{chr}{chromosome name, such as "chr1", "chr22", "chrX" and "chrY", the argument needs to be specified only when bsparse = TRUE.}
\item{resolution}{resolution of Hi-C matrix, the argument needs to be specified only when bsparse = TRUE.}
}
\value{
A dense contact matrix.
}
\description{
Load a Hi-C matrix file in sparse or dense format, and output the dense contact matrix.
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
}
