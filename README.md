### CBLRR

CBLRR: A Cachy-Based Bounded Constraint Low-Rank Representation Method for Single-cell RNA-seq Data Analysis

### Requires: Matlab/octave 6.4.0


### 1. Run CBLRR

To apply CBLRR, please run the Matlab script file `demo.m`. We provide real examples used in our manuscript (Kolod) to reproduce our results. 


### 2. Input Arguments
Please see the /Data/Examples/Input folder and the script file `demo.m` for details.

* `in_X.mat`: the expression matrix of scRNA-seq data used in the CBLRR method. Taking the kolod data matrix as an example, the rows are the genes (10685) and the columns are the samples (704).

* `true_labs.mat`: the labels of scRNA-seq data used in the CBLRR method. Taking the Kolod labels as an example, it contains 5 types (1, 2, ..., 5).

* `n_space`: the number of clusters, predefined by the number of eigenvalues. The default setting is the number of cell types if the number of clusters is given. 

* `alpha`: the tuning parameter for the bounded nuclear norm regulation. We set the value of this parameter by default as 1.

* `a`: the scale parameter for regularizing the cauchy loss function. We set the value of this parameter by default as 10.

* `mu`: the tuning parameter for optimizing and solving the objective function. We set the value of this parameter by default as 110.

* `beta`: the tuning parameter for optimizing and solving the objective function. We set the value of this parameter by default as 1.


### 3. Output variables

There are five output variables. Please see the /Data/Examples/Output folder and the following details.

* `NMI`: the Normalized Mutual Information of clustering results.

* `ARI`: the Adjusted Rand Index of clustering results.

* `grps`: the cell labels predicted by the CBLRR clustering method.

* `Similarity`: the local similarity matrix to learn the heterogeneity of cells and can be used to indicate cell clusters, formed as matrix format.

* `Z`: the coefficient matrix obtained by the CBLRR method, formed as matrix format.


### 4. Files:

CBLRR.m - The main function.

demo.m - A script with a real scRNA-seq data to show how to run the code.

LRR.m - Solve the low-rank problem via ADMM.

Bound.m - Compute the bounded constraint.

Contingency.m - Form contigency matrix for two vectors.

FilterGenesZero.m - Count the zero expression of genes in different cells, and filter genes to control the data quality.

SpectralClustering.m - Computes the clustering of the nodes using the spectral clustering algorithm.

cal_eigenvalues.m - Calculate the number of eigenvalues to determine the number of clusters(If the number of clusters is not given).

data file: Due to the space limitation of github, we only give some datasets, other datasets can be downloaded from (https://hemberg-lab.github.io/scRNA.seq.datasets/) or (https://github.com/10XGenomics/single-cell-3prime-paper) and (https://doi.org/10.6084/m9.figshare.5829687.v7). While Examples file provides the example of input data and output data.

### 5. Example:

Follow the steps below to run CBLRR (also contained in the " demo.m" file). Here use a real scRNA-seq dataset (Kolod) as an example.

clear all;
clc;

%% input:
load('Data\Examples\Input\in_X.mat') %Loading data, each column denotes a gene and each row denotes a cell.
load('Data\Examples\Input\true_labs.mat') %Loading labels.

%% set tuning parameters:
alpha=1; % set as 1 by default.
mu =110; % set as 110 by default.
a = 10; % set as 10 by default.
beta = 1; %set as 1 by default.
n_space = length(unique(true_labs));% The cluster is predefine:
%[n_space,eigenvalues] = cal_eigenvalues(Z);% The number of cluster is not given:

%% perform CBLRR:
[NMI,ARI,grps,similarity,Z] = CBLRR(in_X,true_labs,n_space,alpha,beta,mu,a);


### If you have any questions, please contact the author Qian Ding(dingqian19@126.com).
