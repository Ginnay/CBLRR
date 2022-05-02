clc
clear

load('in_X.mat') %Loading data, each column denotes a gene and each row denotes a cell.
load('true_labs.mat') %Loading labels.

alpha=1;
mu =110;
a = 10;
beta = 1;
n_space = length(unique(true_labs));% The cluster is predefine:

%[n_space,eigenvalues] = cal_eigenvalues(Z);% The number of cluster is not given:
[NMI,ARI,grps,similarity,Z] = CBLRR(in_X,true_labs,n_space,alpha,beta,mu,a);
