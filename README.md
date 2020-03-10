TADBD: a sensitive and fast method for detection of TAD boundaries

Abstract<br>
 　 Topological Associated Domain (TAD) is a self-interacting genomic block. The detection of TAD boundaries on Hi-C contact matrix is one of the most important issues in the analysis of 3D genome architecture at TAD level. Here, we present TADBD, a sensitive and fast computational method for detection of TAD boundaries on Hi-C contact matrix. It implements a Haar-based algorithm considering Haar diagonal template, the acceleration via a compact integrogram, multi-scale aggregation at template size, and statistical filtering. Compared with the other five methods, including HiCDB, IC-Finder, EAST, TopDom and HiCseg, on simulated and experimental data, TADBD shows a better performance in most cases. In addition, a freely available R package is given to facilitate the use of the presented method.

1> Datasets<br>
  Simulated data:<br>
 　   Yu, W., He, B. and Tan, K. (2017) Identifying topologically associating domains and ubdomains by Gaussian Mixture model And 
      Proportion test, Nat Commun, 8, 535.https://bitbucket.org/mforcato/hictoolscompare/src/master/<br> 
  Real data:<br>
 　   Rao, S.S.P., et al. (2014) A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping, Cell, 159,         1665-1680.<br> 
		
2> TADBD<br>		
  Description:<br>
 　   An R package for detection of TAD boundaries<br> 
  Depends:<br>
 　   R (>= 3.5.1)<br>
 Installation:
  install.packages("devtools") # if you have not installed "devtools" package
  devtools::install_github("bioinfo-lab/TADBD")
  Parameter tunning:
  XXXXXXXX

Example:
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
  Output(df_result, species, chr, resolution, outxtfile="./result", bheatmap = T, heatmapfile="./heatmap", hicmat)
