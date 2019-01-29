function [adj_mat, clusters] = standard_SBM(n,r,p,q)

% Dec 2018
% This matlab code implements the standard stochastic block model for
% generating random graphs with planted partitions.
%
%
% n - number of nodes
%
% r - number of clusters
%
% p - probability of existance of an intra-cluster edge
%
% q - probability of existance of an inter-cluster edge
%
%

k=n/r;
clusters=zeros([n 1]);
for i=1:r
    for j=(i-1)*k+1:i*k
        clusters(j)=i;
    end
end
clusters = clusters(randperm(length(clusters)));

adj_mat=eye(n);
for i=1:n
    for j=i+1:n
        prob=rand;
        if clusters(i)==clusters(j)
            if prob <= p
                adj_mat(i,j)=1;
                adj_mat(j,i)=1;
            end
        else
            if prob <= q
                adj_mat(i,j)=1;
                adj_mat(j,i)=1;
            end
        end
    end
end