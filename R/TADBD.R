# @author : Lv Hongqiang and Li Lin
# @brief : TADBD.R is a software code to identify topological domains for given Hi-C contact matrix.
# @version  1.0
# @fn HMDTB
# input
# @param hicmat : matrix, matrixFile by hicdata
# @param template.sizes :vector, some numbers of bins to extend
#                        template.sizes, default value c(4, 5, 6)
#                        template.sizes = others, such as c(3, 4, 5, 6) and so on
# @param bstatfilter : string, bstatfilter = F, not statistical filter
#                              bstatfilter = T, get qvalu and statistical filter (default)
# output
# @param df_result : dataframe, include coordinate information of boundary points
### Supplementary notes ###
# The source of function findpeaks is Rpackage pracma::findpeaks.
# The source of function qvalue is a software package in Bioconductor.
TADBD<-function(hicmat, template.sizes = c(4,5,6), bstatfilter = T)
{
  disdata1 <- hicmat
  size_diag_max <- max(template.sizes)
  distance <- 2*size_diag_max+1

  zeroM <- matrix(0, nrow=nrow(disdata1)+1, ncol=ncol(disdata1)+1)
  zeroM[2:nrow(zeroM), 2:ncol(zeroM)] <- disdata1

  my_matrix_integral <- get_integral(zeroM,distance)

  # Select the window size and do the arithmetic average
  haar_value_diag_add <- getscore(my_matrix_integral, 2)
  haar_value_diag_add <- c(matrix(0, nrow=1, ncol=length(haar_value_diag_add)))

  for (i in 1 : length(template.sizes))
  {
    haar_value_diag_num <- getscore(my_matrix_integral, template.sizes[i])
    haar_value_diag_add <- haar_value_diag_add + haar_value_diag_num
  }
  haar_value_diag <- haar_value_diag_add/length(template.sizes)

  x_peaks<- c(findpeaks(haar_value_diag, nups = 1, minpeakdistance = 2)[,2])
  x_peaks <- sort(x_peaks)
  pvalue<-getpvalue(disdata1,x_peaks,size_diag_max)

  # Caculate qvalue
  qvalue<-qvalue(pvalue)
  qvalue<-qvalue$qvalues

  from.id <- c(0, (x_peaks - 1)[1:(length((x_peaks-1))-1)])
  to.id <- c((x_peaks - 1))

  pvalue<-c(pvalue)
  qvalue<-c(qvalue)


  if (bstatfilter == T)
  {
    from.id1 <- from.id
    to.id1 <- to.id
    outfile1 <- data.frame(from.id1, to.id1, pvalue, qvalue)
    j=1
    boundary<-c()
    for(i in 1:length(x_peaks))
    {
      if(qvalue[i]< 0.05)
      {
        boundary[j]<-x_peaks[i]
        j <- j + 1
      }
    }
    from.id2 <- c(0, (boundary - 1)[1:(length((boundary-1))-1)])
    to.id2 <- c((boundary - 1))
    outfile2 <- data.frame(from.id2, to.id2)
    df_result <- list("allpeaks" = outfile1, "domains" = outfile2)
    return(df_result)
  }
  if (bstatfilter == F)
  {
    from.id3 <- from.id
    to.id3 <- to.id
    outfile3 <- data.frame(from.id3, to.id3)
    df_result <- list("allpeaks" = outfile3)
    return(df_result)
  }
}

############ compact integrogram ##############
get_integral<-function(data,distance)
{

  zeroM<-data
  my_matrix_integral <- matrix(data=0,nrow=nrow(zeroM),ncol=ncol(zeroM))
  for (j in 2:(nrow(my_matrix_integral)-distance+1)) {
    for (i in j:(j+distance-1)) {
      if(i == j)
      {
        my_matrix_integral[i,j] <- 2*my_matrix_integral[i,j-1] - my_matrix_integral[i-1,j-1] + zeroM[i,j]
      }
      else if(i == (j+distance-1))
      {
        my_matrix_integral[i,j] <- my_matrix_integral[i-1,j] + zeroM[i,j]
      }
      else
      {
        my_matrix_integral[i,j] <- my_matrix_integral[i-1,j] - my_matrix_integral[i-1,j-1] + my_matrix_integral[i,j-1] + zeroM[i,j]
      }
      my_matrix_integral[j,i] <- my_matrix_integral[i,j]
    }
  }

  for (j in (nrow(my_matrix_integral)-distance+2):nrow(my_matrix_integral)) {
    for (i in j:nrow(my_matrix_integral)) {
      if(i == j)
      {
        my_matrix_integral[i,j] <- 2*my_matrix_integral[i,j-1] - my_matrix_integral[i-1,j-1] + zeroM[i,j]
      }
      else
      {
        my_matrix_integral[i,j] <- my_matrix_integral[i-1,j] - my_matrix_integral[i-1,j-1] + my_matrix_integral[i,j-1] + zeroM[i,j]
      }
      my_matrix_integral[j,i] <- my_matrix_integral[i,j]
    }
  }
  return(my_matrix_integral)
}

############# get Haar feature value #############
getscore<-function(my_matrix_integral,size_diag)
{
  step <- 1
  i <- 1+size_diag
  j <- 1
  haar_value_diag <- matrix(data=NA,nrow=1)
  while(i <= (nrow(my_matrix_integral)-size_diag))
  {
    a1 <- my_matrix_integral[i,i]
    b1 <- my_matrix_integral[i-size_diag,i-size_diag]
    c1 <- my_matrix_integral[i,i-size_diag]
    d1 <- my_matrix_integral[i-size_diag,i]
    a2 <- my_matrix_integral[i+size_diag,i+size_diag]
    c2 <- my_matrix_integral[i+size_diag,i]
    d2 <- my_matrix_integral[i,i+size_diag]
    c3 <- my_matrix_integral[i+size_diag,i-size_diag]
    d4 <- my_matrix_integral[i-size_diag,i+size_diag]
    value <- (4*a1+b1-2*c1-2*d1+a2-2*c2-2*d2+c3+d4)/(size_diag*size_diag)
    haar_value_diag[j] <- c(value)
    j <- j+1
    i <- i+step
  }
  fill_diag <- matrix(0,ncol = size_diag)

  haar_value_diag <- c(fill_diag,haar_value_diag,fill_diag)
  haar_value_diag <- haar_value_diag/max(haar_value_diag)

  return(haar_value_diag)
}

################# p-value ##################
getpvalue <- function(mat,peaks,size_diag)
{
  pvalue <- rep(1,length(peaks))
  left_up<-matrix(nrow=size_diag+1,ncol=size_diag+1)
  left_down<-matrix(nrow=size_diag+1,ncol=size_diag+1)
  right_up<-matrix(nrow=size_diag+1,ncol=size_diag+1)
  right_down<-matrix(nrow=size_diag+1,ncol=size_diag+1)

  # The data of the sub-diagonal are normalized
  n_bins<-nrow(mat)
  for(i in 0:(2*size_diag+1))
  {
    mat[seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins)] <- scale(mat[seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins)]) #log10(mat[seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins)]+1)
    mat[seq((1+i), n_bins*(n_bins-i), 1+n_bins)] <- mat[seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins)]#scale(mat[ seq((1+i), n_bins*(n_bins-i), 1+n_bins)])
  }
  scale.data<-matrix(0,nrow=(n_bins+2*size_diag),ncol=(n_bins+2*size_diag))
  scale.data[(size_diag+1):(n_bins+size_diag),(size_diag+1):(n_bins+size_diag)] <- mat

  # Do rank sum test
  for(i in 1:length(peaks))
  {
    left <- peaks[i]-size_diag+size_diag
    mid <- peaks[i]+size_diag
    right <- peaks[i]+size_diag-1+size_diag
    left_up<-scale.data[left:(mid-1),left:(mid-1)]
    right_down<-scale.data[mid:right,mid:right]
    right_up<-scale.data[left:(mid-1),mid:right]
    left_down<-scale.data[mid:right,left:(mid-1)]
    left_up<-as.vector(left_up)
    left_down<-as.vector(left_down)
    right_up<-as.vector(right_up)
    right_down<-as.vector(right_down)
    wil.test =  wilcox.test(x=c(right_up,left_down), y=c(left_up,right_down), alternative="less", exact=F)
    pvalue[i] = wil.test$p.value
  }
  return(pvalue)
}

############## q-value ##################
qvalue <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...)
{
  # Argument checks
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  }

  # Calculate pi0 estimate
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1)  {
      pi0s = list()
      pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }

  # Calculate q-value estimates
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 - (1 - p[o]) ^ m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  }
  qvals_out[rm_na] <- qvals
  # Calculate local FDR estimates
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }

  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda,
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}
pi0est <- function(p, lambda = seq(0.05,0.95,0.05), pi0.method = c("smoother", "bootstrap"),
                   smooth.df = 3, smooth.log.pi0 = FALSE, ...)
{
  # Check input arguments
  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method = match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda) # guard against user input

  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  } else if (ll > 1 && ll < 4) {
    stop(sprintf(paste("ERROR:", paste("length(lambda)=", ll, ".", sep=""),
                       "If length of lambda greater than 1,",
                       "you need at least 4 values.")))
  } else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }

  if (max(p) < max(lambda)) {
    stop("ERROR: maximum p-value is smaller than lambda range. Change the range of lambda or use qvalue_truncp() for truncated p-values.")
  }

  # Determines pi0
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  } else {
    ind <- length(lambda):1
    pi0 <- cumsum(tabulate(findInterval(p, vec=lambda))[ind]) / (length(p) * (1-lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    # Smoother method approximation
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      } else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    } else if (pi0.method == "bootstrap") {
      # Bootstrap method closed form solution by David Robinson
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W / (m ^ 2 * (1 - lambda) ^ 2)) * (1 - W / m) + (pi0 - minpi0) ^ 2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    } else {
      stop('ERROR: pi0.method must be one of "smoother" or "bootstrap".')
    }
  }
  if (pi0 <= 0) {
    warning("The estimated pi0 <= 0. Setting the pi0 estimate to be 1. Check that you have valid p-values or use a different range of lambda.")
    pi0 <- pi0.lambda <- 1
    pi0Smooth <- lambda <- 0
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda,
              lambda = lambda, pi0.smooth = pi0Smooth))
}
lfdr <- function(p, pi0 = NULL, trunc = TRUE, monotone = TRUE,
                 transf = c("probit", "logit"), adj = 1.5, eps = 10 ^ -8, ...)
{
  # Check inputs
  lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range [0,1].")
  } else if (is.null(pi0)) {
    pi0 <- pi0est(p, ...)$pi0
  }
  n <- length(p)
  transf <- match.arg(transf)
  # Local FDR method for both probit and logit transformations
  if (transf == "probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)
    x <- qnorm(p)
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0 * dnorm(x) / y
  } else {
    x <- log((p + eps) / (1 - p + eps))
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1 + exp(x)) ^ 2
    lfdr <- (pi0 * dx) / y
  }
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {
    o <- order(p, decreasing = FALSE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  return(lfdr_out)
}
########### find peaks #################
findpeaks <- function(x,nups = 1, ndowns = nups, zero = "0", peakpat = NULL,
                      # peakpat = "[+]{2,}[0]*[-]{2,}",
                      minpeakheight = -Inf, minpeakdistance = 1,
                      threshold = 0, npeaks = 0, sortstr = FALSE)
{
  stopifnot(is.vector(x, mode="numeric") || length(is.na(x)) == 0)
  if (! zero %in% c('0', '+', '-'))
    stop("Argument 'zero' can only be '0', '+', or '-'.")

  # transform x into a "+-+...-+-" character string
  xc <- paste(as.character(sign(diff(x))), collapse="")
  xc <- gsub("1", "+", gsub("-1", "-", xc))
  # transform '0' to zero
  if (zero != '0') xc <- gsub("0", zero, xc)

  # generate the peak pattern with no of ups and downs
  if (is.null(peakpat)) {
    peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
  }

  # generate and apply the peak pattern
  rc <- gregexpr(peakpat, xc)[[1]]
  if (rc[1] < 0) return(NULL)

  # get indices from regular expression parser
  x1 <- rc
  x2 <- rc + attr(rc, "match.length")
  attributes(x1) <- NULL
  attributes(x2) <- NULL

  # find index positions and maximum values
  n <- length(x1)
  xv <- xp <- numeric(n)
  for (i in 1:n) {
    xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
    xv[i] <- x[xp[i]]
  }

  # eliminate peaks that are too low
  inds <- which(xv >= minpeakheight & xv - pmax(x[x1], x[x2]) >= threshold)

  # combine into a matrix format
  X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])

  # eliminate peaks that are near by
  if (minpeakdistance < 1)
    warning("Handling 'minpeakdistance < 1' is logically not possible.")

  # sort according to peak height
  if (sortstr || minpeakdistance > 1) {
    sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
    X <- X[sl, , drop = FALSE]
  }

  # return NULL if no peaks
  if (length(X) == 0) return(c())

  # find peaks sufficiently distant
  if (minpeakdistance > 1) {
    no_peaks <- nrow(X)
    badpeaks <- rep(FALSE, no_peaks)

    # eliminate peaks that are close to bigger peaks
    for (i in 1:no_peaks) {
      ipos <- X[i, 2]
      if (!badpeaks[i]) {
        dpos <- abs(ipos - X[, 2])
        badpeaks <- badpeaks | (dpos > 0 & dpos < minpeakdistance)
      }
    }
    # select the good peaks
    X <- X[!badpeaks, , drop = FALSE]
  }

  # Return only the first 'npeaks' peaks
  if (npeaks > 0 && npeaks < nrow(X)) {
    X <- X[1:npeaks, , drop = FALSE]
  }

  return(X)
}
