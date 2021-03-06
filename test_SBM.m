%
% Dec 2018
% This matlab code is to test the Improved Graph Clustering method, on
% graphs generated by standard SBM
%
%
%

addpath SRC;

n=1000;
% number of clusters should be fixed for fairness of comparison.
r=5;

ALM_p_avg=zeros([1 13]);
slink_p_avg=zeros([1 13]);
theory_p=zeros([1 13]);
step1=1.05;
step2=1.1;
i=0;
for q=0:0.05:0.6
    syms p positive
    eqn = p-q-sqrt(p*(1-q)*n)/(n/r);
    solp = vpasolve(eqn,p);
    p = double(solp)+0.02;
    i=i+1;
    disp(['#round ' num2str(i) ': q=' num2str(q) ' theory_p=' num2str(p)])
    for j=1:20
        ALM_pur=0;
        ALM_p=p;
        while ALM_pur<0.99 && ALM_p<1
            ALM_p=step1*ALM_p;
            [A, clusters] = standard_SBM(n,r,ALM_p,q);
            [ALM_cluster, A_dual] = improved_graph_cluster(A,r);
            ALM_pur=purity(clusters, ALM_cluster);
        end
        if ALM_p>1
            ALM_p=1;
        end
        ALM_p_avg(i) = ALM_p_avg(i) + ALM_p;
        
        slink_pur=0;
        slink_p=p;
        while slink_pur<0.99 && slink_p<1
            slink_p=step2*slink_p;
            [A, clusters] = standard_SBM(n,r,slink_p,q);
            slink_tree = linkage(A, 'single');
            slink_cluster = cluster(slink_tree,'Maxclust',r);
            slink_pur=purity(clusters, slink_cluster);
        end
        if slink_p>1
            slink_p=1;
        end
        slink_p_avg(i) = slink_p_avg(i) + slink_p;
        
        disp(['    #trial ' num2str(j) ': IGC_p=' num2str(ALM_p)...
            ' SLINK_p=' num2str(slink_p)])
    end
    ALM_p_avg(i)=ALM_p_avg(i)/j;
    slink_p_avg(i) = slink_p_avg(i)/j;
    theory_p(i)=p;
end

figure()
plot(0:0.05:0.6,slink_p_avg,'o-')
hold on
plot(0:0.05:0.6,ALM_p_avg,'^-')
plot(0:0.05:0.6,theory_p,'.-')
hold off
legend('SLINK','IGC','Theory')
title('avg minimum p satisfying purity > 0.99 over 20 trials')
xlabel('q')
ylabel('p')