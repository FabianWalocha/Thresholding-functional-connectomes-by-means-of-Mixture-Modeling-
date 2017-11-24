Please read carefully before using the methods in this repository. This set of codes was intended to implement the methods depicted in the paper Thresholding functional connectomes by means of Mixture Modeling (2017) by Bielczyk, Walocha et al. 

The main method Mixture Modeling is realized via the function em_alg on a single subject level, whereas the functions do_MM and do_MM2 implement the data on a set of BOLD time series data and a set of partial correlation values respectively. Any method labeled plotXXX is used to create the plots from the paper mentioned above.

All code, unless explicitly stated otherwise, was written by Fabian Walocha (fawalocha@gmail.com) in 2016-2017.

----------------
do_MM - 
	Calculates the sparsified partial correlation matrix using a
	mixture modeling based approach from BOLD time series input datasets. 	
	If no mmtype and lwFlag is given, optimal model is chosen based on 
	Bayesian Information Criterion (BIC) 

do_MM2 - 
	Employs mixture modeling based approach on a set of partial
	correlation matrices
	If no mmtype is given, optimal model is chosen based on Bayesian
	Information Criterion (BIC) 
----------------

em_alg - 
	Does mixture fit using EM algorithm

do_BIC - 
	Calculates the Bayesian Information Criterion for a given set
	of partial correlation matrices and a mixture type

covCor - 
	All credits go to Copyright (c) 2014, Olivier Ledoit and Michael Wolf 

ICC_new - 
	All credits go to Kevin Brownhill, Imaging Sciences, KCL, London kevin.brownhill@kcl.ac.uk
	Taken from mathworks.com/matlabcentral/fileexchange/21501-intraclass-correlation-coefficients?focused=5143573&tab=function on 2017-08-04

saveTightFigure - 
	All credits go to E Akbas (c) Aug 2010
	Taken from mathworks.com/matlabcentral/fileexchange/48877-savetightfigure-h-outfilename- on 2017-09-16

plotBIC - 
	Produces: Figure 6, figure 7, figure 8

plotFDR - 
	Produces: Figure 1C

plotMain0 - 
	Arranges data in order to prepare for plotMain1

plotMain1 - 
	Produces: Figure 1A, figure 1B, Figure 3, table 1, table 2, table 3,
	table 4, saves results for figure 9, figure 10

plotMain2 - 
	Saves codes for figure 4, figure 5, figure 11

plotROC - 
	Produces: Figure 1D

plotVarN - 
	Produces: Figure 2

circularPlots - 
	Produces: Figure 4, figure 5, figure 9, figure 10, figure 11
	Made using python 3.5.3.1Qt5