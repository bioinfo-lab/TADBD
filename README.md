TADBD: a sensitive and fast method for detection of TAD boundaries

Abstract<br>
 　  Topological Associated Domain (TAD) is a self-interacting genomic block. The detection of TAD boundaries on Hi-C contact matrix is one of the most important issues in the analysis of 3D genome architecture at TAD level. Here, we present TADBD, a sensitive and fast computational method for detection of TAD boundaries on Hi-C contact matrix. It implements a Haar-based algorithm considering Haar diagonal template, the acceleration via a compact integrogram, multi-scale aggregation at template size, and statistical filtering. Compared with the other six popular methods based on simulated and experimental Hi-C data, TADBD is competitive in terms of accuracy, speed and robustness. In addition, a new R package TADBD is freely given out online.

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
  Example:
 
    rm(list=ls())
	library(TADBD)
	species <- "hg19"
	chr <- "chr18"
	resolution <- 50000
	options(scipen = 999)
	data(hicdata)
	hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)
	df_result <- TADBD(hicmat)
	#Output a text file
	Output(df_result, species, chr, resolution)
<<<<<<< HEAD
	#Output a heatmap
	Output(df_result, species, chr, resolution, outxtfile="./result", bheatmap = T, heatmapfile="./heatmap", hicmat)
=======
	//: Output a heatmap
    Output(df_result, species, chr, resolution, outxtfile="./result", bheatmap = T, heatmapfile="./heatmap", hicmat)  

>>>>>>> ce5578511a2204f33303bc211d58b5266fa38209
         
