\name{TADBD}
\alias{TADBD}
\title{Take a dense contact matrix as input, and output the bin number of TAD
       boundaries on the contact matrix}
\usage{
TADBD(hicmat, template.sizes, bstatfilter)

hicmat            a contact matrix by DataLoad function.
template.sizes    a vector consisting of multiple template sizes, with a default
                  value of c(4, 5, 6).The vector can also be specified by user, for
                  example: c(3, 4, 5, 6).
bstatfilter       if statistical filtering is involved, bstatfilter = T (default),
                  otherwise bstatfilter = F.
}
\value{
A dataframe consisting of the bin number of the detected TAD boundaries.
}
\description{
Take a dense contact matrix as input, and output the bin number of TAD boundaries on the contact matrix.
}
\examples{
#Load R package TADBD
library(TADBD)
#Configuration of the parameters, including species, chromsome and resolution
species <- "hg19"
chr <- "chr18"
resolution <- 50000
#Specify Hi-C data to be loaded
options(scipen = 999)
#Specify Hi-C data to be loaded
data(hicdata)
#Load a Hi-C contact matrix file in a dense format
hicmat <- DataLoad(hicdata, bsparse = FALSE, species, chr, resolution)
#Detect TAD boundaries on the loaded contact matrix using TADBD method
df_result <- TADBD(hicmat)
}
