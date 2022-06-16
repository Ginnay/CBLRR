clc
clear

%% input:
% The expression matrix of scRNA-seq data used in the CBLRR method.Taking the kolod data matrix as an example, the rows are the genes (10685) and the columns are the samples (704).
% The labels of scRNA-seq data used in the CBLRR method. Taking the Kolod labels as an example, it contains 5 types (1, 2, ..., 5).
load('\Data\Examples\Input\in_X.mat') %Loading data, each column denotes a gene and each row denotes a cell.
load('\Data\Examples\Input\true_labs.mat') %Loading labels.

%% set tuning parameters:
alpha=1; % set as 1 by default.
mu =110; % set as 110 by default.
a = 10; % set as 10 by default.
beta = 1; %set as 1 by default.
n_space = length(unique(true_labs));% The cluster is predefine:
%[n_space,eigenvalues] = cal_eigenvalues(Z);% The number of cluster is not given:

%% perform CBLRR:
[NMI,ARI,grps,similarity,Z] = CBLRR(in_X,true_labs,n_space,alpha,beta,mu,a);

