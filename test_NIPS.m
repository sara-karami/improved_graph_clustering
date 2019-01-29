%
% Dec 2018
% This matlab code is to test the Improved Graph Clustering method, on NIPS
% dataset, available at: https://cs.nyu.edu/~roweis/data.html
%
%
%

addpath SRC;

A=load('nipsAdjMatrix.csv');
% number of clusters is fixed to r=8 for fairness of comparison.
ALM_in_density = zeros([1 5]);
ALM_cross_density = zeros([1 5]);
slink_in_density = zeros([1 5]);
slink_cross_density = zeros([1 5]);
i=1;
for r=6:10
ALM_cluster = improved_graph_cluster(A,r);
[ALM_in_density(i), ALM_cross_density(i)] = cluster_density(ALM_cluster, A, r);

slink_tree = linkage(A, 'single');
slink_cluster = cluster(slink_tree,'Maxclust',r);
[slink_in_density(i), slink_cross_density(i)] = cluster_density(slink_cluster, A, r);

disp(['#clusters = ' num2str(r)]);
disp(['    IGC: intra_cluster_density=' num2str(ALM_in_density(i))...
    ' inter_cluster_density=' num2str( ALM_cross_density(i))])
disp(['    SLINK: intra_cluster_density=' num2str(slink_in_density(i))...
    ' inter_cluster_density=' num2str( slink_cross_density(i))])
i=i+1;
end

figure()
subplot(1,2,1)
plot(6:r,slink_in_density,'o-')
hold on
plot(6:r,ALM_in_density,'^-')
hold off
legend('SLINK','IGC')
xlabel('r')
ylabel('intra cluster density')
grid on
subplot(1,2,2)
plot(6:r,slink_cross_density,'o-')
hold on
plot(6:r,ALM_cross_density,'^-')
hold off
legend('SLINK','IGC')
xlabel('r')
ylabel('inter cluster density')
grid on