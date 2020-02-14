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
rm(list=ls())
library(TADBD)
species <- "hg19"
chr <- "chr18"
resolution <- 50000
options(scipen = 999)
data(hicdata)
hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)
df_result <- TADBD(hicmat)
}
