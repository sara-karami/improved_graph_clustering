function [r,k,p,q,t] = estimate_param(A,fixed_r)

% Dec 2018
% This matlab code computes clustering parameters needed for Improved Graph
% clustering method.
%
% C - n x 1 vector of node clusters
%
% A - n x n adjacency matrix of input graph
%
% K - number of clusters
%
%
%

e = eig(A);
e = sort(e,'descend');
n = length(A);
if nargin>1 && fixed_r>0
    r=fixed_r;
else
    r = e(2)-e(3);
    for i=2:n-1
        if e(i)-e(i+1)>r
            r = e(i)-e(i+1);
        end
    end
end
r = round(r,0);
% k = cluster size
k = n/r;
% p = probability of within cluster edges
p = (k*e(1)+(n-k)*e(2)-n)/(n*(k-1));
% q = probability of between cluster edges
q = (e(1)-e(2))/n;
t = (p+q)/2;
end