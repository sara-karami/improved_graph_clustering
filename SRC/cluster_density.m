function [in_cluster_density, cross_cluster_density] = cluster_density(C, A, K)
% Dec 2018
% This matlab code computes in-cluster and cross-cluster edge densities
% These values are used to evaluate clustering quality
%
%
% C - n x 1 vector of node clusters
%
% A - n x n adjacency matrix of input graph
%
% K - number of clusters
%
%
%

n=length(C);
in_cluster_density=0;
max_in_cluster_density=0;
cross_cluster_edge=0;
for c_num=1:K
    c_size=0;
    in_cluster_edge=0;
    for i=1:n
        if C(i)==c_num
            c_size=c_size+1;
        end
        for j=1:n
            if C(i)==c_num && C(j)==c_num
                if A(i,j)==1
                    in_cluster_edge=in_cluster_edge+1;
                end
            elseif C(i)~=C(j)
                if A(i,j)==1
                    cross_cluster_edge=cross_cluster_edge+1;
                end
            end
        end
    end
    max_in_cluster_density=max_in_cluster_density+c_size*(c_size-1);
    if c_size>1
        in_cluster_density=(in_cluster_edge/(c_size*(c_size-1)))/K+in_cluster_density;
    end
    %c_size
end
cross_cluster_density=cross_cluster_edge/(n*(n-1)-max_in_cluster_density);
