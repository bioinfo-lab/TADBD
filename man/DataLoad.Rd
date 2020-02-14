\name{DataLoad}
\alias{DataLoad}
\title{Load a Hi-C matrix file in sparse or dense format, and output the dense
       contact matrix}
\usage{
DataLoad(hicdata, bsparse, species, chr, resolution)

hicdata        a Hi-C matrix file.
bsparse        if the Hi-C matrix is in sparse format, bsparse = T, otherwise
               bsparse = F.
species        species name, including "hg19", "hg38", "mm9" and "mm10", the
               argument needs to be specified only when bsparse = T.
chr            chromosome name, such as "chr1", "chr22", "chrX" and "chrY", the
               argument needs to be specified only when bsparse = T.
resolution     resolution of Hi-C matrix, the argument needs to be specified only
               when bsparse = T.
}
\value{
A dense contact matrix.
}
\description{
Load a Hi-C matrix file in sparse or dense format, and output the dense contact matrix.
}
\examples{
rm(list=ls())
library(TADBD)
species <- "hg19"
chr <- "chr18"
resolution <- 50000
options(scipen = 999)
data(hicdata)
hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)
}
