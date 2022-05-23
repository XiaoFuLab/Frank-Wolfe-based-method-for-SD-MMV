function [theta_hat, Tracking] = runByMERIT(X, K, varargin) 
% ----------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('maxIters', 1000);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    ops = p.Results;

    options = struct;
    options.maxIters = 120;
    options.verbose = 1;
    options.K = K;
    options.estimatedK = 4*K;
    options.innerVerboseInterval = int64(10000);
    options.backend = 'mex';
    options.epsilon = 1e-5;
    options.numThreads = 4;
    options.dualGapThreshold = 0.5; 

    originalX = X;
    [originalV, ~] = eigs(X, K);

    [V, D] = eigs(X, K);
    options.firstIterMethod = 'NONE';

    V = truncateData(V, 0.99);
    X = V';
    options.lambda =  1e-6;

    fprintf('Running MERIT ...\n')
    [hatC, Tracking] = MERIT(X, options);
    Tracking.fwX = X;
    fprintf('Done \n')
    [~, lambdaHat] = maxk(vecnorm(hatC, Inf, 2), K, 1);
    
    W_hat = X(:, lambdaHat);
    H_hat = originalV / W_hat';
    tre = 1e-12;
    H_hat(H_hat<tre) = 0;
    H_hat = normalize_row_l1_(H_hat);
    theta_hat = H_hat;
end

function [ theta_l1 ] = normalize_row_l1_( theta )

d_theta = sum(theta,2);
S = find(d_theta>0);
theta_l1 = theta;
theta_l1(S,:) = bsxfun(@times, theta(S,:), 1./(sum(theta(S,:), 2)));

S1 = find(d_theta==0);
theta_l1(S1, :) = ones(numel(S1), size(theta, 2)) .* mean(theta_l1(S, :));

end
