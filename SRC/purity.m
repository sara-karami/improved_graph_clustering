function pur = purity(true_cluster, test_cluster)

% Dec 2018
% This matlab code is to calculate purity factor of a clustering result
% given the true clusters and the results.
%
%
% true_cluster - n x 1 vector of nodes' true clusters
%
% test_cluster - n x 1 vector of nodes' clustering results
%
% As the order of given clusters has nothing to do with the results so it
% can also be used in this form: purity(cluster1, cluster2) where one of
% the inputs is the true cluster vector and the other one is the test
% cluster vector.
%
%

cross_freq=crosstab(true_cluster,test_cluster);
sum_freq=sum(cross_freq,1);
max_freq=max(cross_freq,[],1);
pur=sum(max_freq(sum_freq>0)/sum(sum_freq));

end