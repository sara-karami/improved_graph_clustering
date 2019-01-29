function [Y_hat, S_hat, iter] = inexact_alm_revised(A, C, lambda, epsilon, max_iter, alpha)

% Dec 2018
% This matlab code implements the weighted inexact augmented Lagrange 
% multipier method for Graph Clustering.
% This code is obtained from altering the ALM code provided by: 
% Zhouchen Lin, Minming Chen, and Yi Ma, The Augmented Lagrange Multiplier Method 
% for Exact Recovery of Corrupted Low-Rank Matrix, http://perception.csl.illinois.edu/matrix-rank/Files/Lin09-MP.pdf
%
%
% A - n x n adjacency matrix (required input)
%
% C - n x n Weight matrix
%
% lambda - weight on sparse error term in the cost function
%
% epsilon - tolerance for stopping criterion
%
% max_iter - maximum number of iterations
%
% mu - learning rate
% 
% alpha - increasing ratio of learning rate
%
% Initialize Y,S,M,mu
% while ~converged 
%   minimize (inexactly, update Y and S only once)
%     L(Y,S,M,mu) = |Y|_* + lambda * |C .* S|_1 + <M,A-Y-S> + mu/2 * |A-Y-S|_F^2;
%   M = M + \mu * (A -Y -S);
%   update \mu; %The next \mu is predicted and clipped to [\mu,\alpha*\mu]
% end
%
%
%

addpath SRC\PROPACK;

[n, ~] = size(A);


%% initialize
M = zeros(n);
Y_hat = zeros(n);
S_hat = zeros(n);
mu = 1.25/norm(A,2);
mu_bar = mu * 1e9;
tol2 = 1e-5;
norm_A = norm(A, 'fro');

iter = 0;
converged = false;
sv = min([10, n]);

%% Weighted IALM
while ~converged
    iter = iter + 1;
    

    %% update Y
    if choosvd(n, sv) == 1
        [U, S, V] = lansvd(A - S_hat + M/mu, sv, 'L');
    else
        [U, S, V] = svd(A - S_hat + M/mu, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min([svp + 1, n]);
    else
        sv = min([svp + round(0.05*n), n]);
    end
    temp1 = Y_hat;
    Y_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';
    
    for i1=1:n
        for j1=1:n
            Y_hat(i1,j1)=max(min(Y_hat(i1,j1),1),0);
        end
    end

    %% update S
    temp = A - Y_hat + M/mu;

for i1=1:n
    for j1=1:n
        S_hat(i1,j1) = max(temp(i1,j1) - lambda/mu*C(i1,j1), 0);
        S_hat(i1,j1) = S_hat(i1,j1) + min(temp(i1,j1) + lambda/mu*C(i1,j1), 0);
    end
end

    %% one stop criterion is mu_k*||Y_{k+1}-Y_k||/||A|| < tol2
    converged2 = false;
    stopCriterion2 = mu*norm(Y_hat - temp1, 'fro') / norm_A;
    if stopCriterion2 < tol2
        converged2 = true;
    end      
    
    %% the other stop criterion is ||A-Y-S||/||A|| < tol  
    converged = false;
    temp = A - Y_hat - S_hat;      
    stopCriterion1 = norm(temp, 'fro') / norm_A;
    if stopCriterion1 < epsilon
        converged = true;
    end    
    
    %% update mu
    if converged2
       mu = min(mu*alpha, mu_bar);
    end
    
    %% update M
    M = M + mu*temp;    
       
    norm_M = mu*norm(M(:),inf)/lambda;
%     if mod(iter, 1) == 0
%         disp(['#svd ' num2str(iter) ' r(A) ' num2str(svp)...
%             ' |E|_0 ' num2str(length(find(abs(S_hat)>0)))...
%             ' stopCriterion1 ' num2str(stopCriterion1)...
%             ' stopCriterion2 ' num2str(stopCriterion2)...
%             ' ||Y|| ' num2str(norm_M)...
%             ' mu ' num2str(mu)]);
%     end    
    
    if ~converged && iter >= max_iter
        %disp('Maximum iterations reached') ;
        converged = 1;              
    end
end
