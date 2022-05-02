#CBLRR

CBLRR: A Cachy-Based Bounded Constraint Low-Rank Representation Method for Single-cell RNA-seq Data Analysis

#Requires: Matlab/octave 6.4.0

#Files:

CBLRR.m - The main function.

demo.m - A script with a real scRNA-seq data to show how to run the code.

in_X.mat -  The expression matrix of a real scRNA-seq data (Kolod) used in the cell type clustering example. The dataset contains 10685 genes and 704 cells.

true_labs.mat - The labels of a real scRNA-seq data (Kolod) used in the cell type clustering example. 

LRR.m - Solve the low-rank problem via ADMM.

Bound.m - Compute the bounded constraint.

Contingency.m - Form contigency matrix for two vectors.

FilterGenesZero.m - Count the zero expression of genes in different cells, and filter genes to control the data quality.

SpectralClustering.m - Computes the clustering of the nodes using the spectral clustering algorithm.

cal_eigenvalues.m - Calculate the number of eigenvalues to determine the number of clusters(If the number of clusters is not given).

data file: Due to the space limitation of github, we only give some datasets, other datasets can be downloaded from (https://hemberg-lab.github.io/scRNA.seq.datasets/) or (https://github.com/10XGenomics/single-cell-3prime-paper) and (https://doi.org/10.6084/m9.figshare.5829687.v7).

#Example:

Follow the steps below to run CBLRR (also contained in the " demo.m" file). Here use a real scRNA-seq dataset (Kolod) as an example.

clear all;

clc;

load('in_X.mat') %Loading data, each column denotes a gene and each row denotes a cell.

load('true_labs.mat') %Loading labels.

alpha=1;
mu =110;
a = 10;
beta = 1;
n_space = length(unique(true_labs));% The cluster is predefine.

%[n_space,eigenvalues] = cal_eigenvalues(Z);% The number of cluster is not given.

[NMIM,ARIM,NMI,ARI,grps,similarity,Z] = CBLRR(in_X,true_labs,n_space,alpha,beta,mu,a);

# If you have any questions, please contact the author Qian Ding(dingqian19@126.com).


