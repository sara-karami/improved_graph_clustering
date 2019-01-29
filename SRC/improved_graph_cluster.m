function [clusters, A_dual] = improved_graph_cluster(A,fixed_r)

% Dec 2018
% This matlab code implements the Improved Graph Clustering method,
% provided by: 
% Yudong Chen, Sujay Sanghavi, and Huan Xu, Improved Graph Clustering,
% https://doi.org/10.1109/TIT.2014.2346205
%
%
% A - n x n adjacency matrix (required input)
%
% fixed_r - number of clusters if it is preferred to be a fixed value
%
%
%


A_mean = mean(mean(A));
[r,k,p,q,t]=estimate_param(A,fixed_r);
n = length(A);

c_A = ((1-t)/t)^0.5; c_Ac = (t/(1-t))^0.5;
C = zeros(n);
for i=1:n
    for j=1:n
        if A(i,j)>A_mean
            C(i,j)=c_A;
        else
            C(i,j)=c_Ac;
        end
    end
end

alpha = 1.5;
lambda = 1/sqrt(n); % it didn't work with 1/(48*sqrt(n))
epsilon = 0.01;
max_iter = 1000;

[A_dual, E_dual, numIter] = inexact_alm_revised(A, C, lambda, epsilon, max_iter, alpha) ;

A_mean=mean(mean(A_dual));
for i=1:n
    for j=1:n
        if A_dual(i,j)>A_mean
            A_dual(i,j)=1;
        else
            A_dual(i,j)=0;
        end
    end
end


ALM_tree = linkage(A_dual,'single');
clusters = cluster(ALM_tree,'Maxclust',r);
end