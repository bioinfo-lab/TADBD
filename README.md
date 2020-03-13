TADBD: a sensitive and fast method for detection of TAD boundaries

Abstract<br>
 　 Topological Associated Domain (TAD) is a self-interacting genomic block. The detection of TAD boundaries on Hi-C contact matrix is one of the most important issues in the analysis of 3D genome architecture at TAD level. Here, we present TADBD, a sensitive and fast computational method for detection of TAD boundaries on Hi-C contact matrix. It implements a Haar-based algorithm considering Haar diagonal template, the acceleration via a compact integrogram, multi-scale aggregation at template size, and statistical filtering. Compared with the other five methods, including HiCDB, IC-Finder, EAST, TopDom and HiCseg, on simulated and experimental data, TADBD shows a better performance in most cases. In addition, a freely available R package is given to facilitate the use of the presented method.

1> Datasets<br>
  Simulated data:<br>
 　   Yu, W., He, B. and Tan, K. (2017) Identifying topologically associating domains and ubdomains by Gaussian Mixture model And 
      Proportion test, Nat Commun, 8, 535.https://bitbucket.org/mforcato/hictoolscompare/src/master/<br> 
  Real data:<br>
 　   Rao, S.S.P., et al. (2014) A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping, Cell, 159, 1665-1680.<br> 
		
2> TADBD<br>		
  Description:<br>
 　   An R package for detection of TAD boundaries<br> 
  Depends:<br>
      Windows 10/linux<br>
      R (>= 3.6.2)<br>
     
 Installation:<br>
 ```
  install.packages("devtools") # if you have not installed "devtools" package
  devtools::install_github("bioinfo-lab/TADBD")
 ```
  Parameter tunning:<br>
  template.sizes: a vector consisting of multiple template sizes with a default value of c(4, 5, 6). It is worth noting that as template sizes increase, the average size of the detected TADs increases and the number of TADs decreases. The template sizes are recommended to be close to (or slightly larger than) the minimal size of expected TADs to achieve higher accuracy.<br>
  bstatfilter: an indicative parameter indicating whether to call statistical filtering, which is usually turned on by default<br>

Example:<br>
```
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
  hicmat <- DataLoad(hicdata, bsparse = FASLE, species, chr, resolution)
  #Detect TAD boundaries on the loaded contact matrix using TADBD method using default parameter configuration, that is template.sizes = c(4,5,6), bstatfilter = TRUE
  df_result <- TADBD(hicmat)
  #Output two text files, one is for detected TAD boundaries, the other for intermediate peaks
  Output(df_result, species, chr, resolution, outxtfile="./result")
  #Output two text files and a heatmap with TAD boundary tracks, the parameters of heatmap include starting and ending coordinates, as well as the color and the width of tracks
  Output(df_result, species, chr, resolution, outxtfile="./result", bheatmap = TRUE, heatmapfile="./heatmap", hicmat, map_start=0, map_end=10000000, l_color="bule", l_width=2.5)
```
