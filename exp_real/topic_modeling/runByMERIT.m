function [K_hat, C_hat, Tracking] = runByMERIT(X, K, varargin) 
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('numOutliers', 0);

    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    ops = struct;
    ops.maxIters = 100;
    ops.verbose = true;
    ops.lambda = 1e-6; 
    ops.K = K;
    ops.estimatedN = 4*K;
    ops.backend = 'mex';
    ops.epsilon = 1e-5;
    ops.dualGapThreshold = 0.5;
    ops.innerVerboseInterval = 10000;

    N = size(X, 2);
    K_hat_spa = SPAselect(X, K);
    W_hat_spa = X(:, K_hat_spa);
    initC = sparse(N, N);
    H_hat_spa = constrained_ls(W_hat_spa, X, 'max_iters', 100);
    initC(K_hat_spa, :) = H_hat_spa;
    ops.initC = initC;
    ops.iterStartAt = 1/sqrt(1/N * norm(X - W_hat_spa*H_hat_spa)^2);
    [C_hat, Tracking] = MERIT(X, ops);
    [~, K_hat] = maxk(vecnorm(C_hat, Inf, 2), N, 1);
end
