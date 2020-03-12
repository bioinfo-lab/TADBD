# @fn DataLoad
# input
# @param hicdata : string, the path of input data
# @param bsparse : string, if input data is density, sparse = FALSE
#                          if input data is sparse,  sparse = TRUE
# @param species : string, species name
# @param chr : string, chromosome number of input data
# @param resolution : int, resolution of input data
# output
# @param hicmat : matrix, matrixFile by hicdata
DataLoad <- function(hicdata, bsparse, species, chr, resolution)
{
  if (bsparse == FALSE)
  {
    #######  not sparse  ####################
    # hicdata<-read.table(hicdata_path, head=FALSE)
    hicmat<-as.matrix(hicdata)
    return(hicmat)
  }
  if (bsparse == TRUE)
  {
    #########  sparse  ######################
    switch(species,
           hg19=
             {
               chrlens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                           146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                           102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                           155270560, 59373566)
               names(chrlens) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                                   "chrX", "chrY")
               dim <- floor(chrlens[chr] / resolution) + 1
             },
           hg38=
             {
               chrlens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
                           145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
                           101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
                           156040895, 57227415)
               names(chrlens) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                                   "chrX", "chrY")
               dim <- floor(chrlens[chr] / resolution) + 1
             },
           mm9=
             {
               chrlens = c(197195432, 181748087, 159599783, 155630120, 152537259, 149517037,	152524553,
                           131738871,	124076172,	129993255,	121843856, 121257530,	120284312,	125194864,
                           103494974,	98319150,	95272651,	90772031,	61342430,	166650296,	15902555,	16300)
               names(chrlens) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                                   "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM")
               dim <- floor(chrlens[chr] / resolution) + 1
             },
           mm10=
             {
               chrlens = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459,
                           129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244,
                           104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698, 16299)
               names(chrlens) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                                   "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM")
               dim <- floor(chrlens[chr] / resolution) + 1
             }
    )

    # hictable<-read.table(hicdata_path, head=FALSE)
    hicmat <- sparse2matrix(hicdata, dim, resolution)
    return(hicmat)
  }
}

########## sparse to matrix ############
sparse2matrix <- function(sparse, dim, resolution)
{
  sparsetriangular <- sparse;
  if (ncol(sparsetriangular) >= 3)
  {
    bins <- unique(c(sparsetriangular[, 1], sparsetriangular[, 2]))
    bins <- as.numeric(bins)
    bins <- bins[order(bins)]
    bins <- unique(bins)
    bin.size <- min(diff(bins))
    if(bin.size != resolution)
    {
      stop("Resolution assigned by user is not compatible with the corrodinates of bins")
    }
    hicmatrix <- matrix(0, dim, dim)
    rownum <- sparsetriangular[, 1] / bin.size + 1
    colnum <- sparsetriangular[, 2] / bin.size + 1
    hicmatrix[cbind(rownum, colnum)] <- as.numeric(c(sparsetriangular[, 3]))
    hicmatrix[cbind(colnum, rownum)] <- as.numeric(c(sparsetriangular[, 3]))
    hicmatrix[is.na(hicmatrix)] <- 0

    return(hicmatrix)
  }
  else
  {
    stop("The number of columns of a sparse matrix should be not less than 3")
  }
}
