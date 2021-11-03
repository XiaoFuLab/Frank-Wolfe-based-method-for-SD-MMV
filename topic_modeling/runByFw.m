function [lambdaHat, hatC, Tracking] = runByFw(X, N, varargin) 
% ----------------------------------------------------------------------
% 
% Same as fw_select_5, but use l2 norm
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('numOutliers', 0);

    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    % FRANK-WOLFE
    ops = struct;
    ops.maxIters = 100;
    ops.verbose = true;
    
    ops.lambda = 1e-6; 

    ops.N = N;
    ops.objTol = 0.00;
    ops.estimatedN = 4*N;
    ops.innerVerboseInterval = int64(5000);
    ops.backend = 'mex';
    ops.epsilon = 1e-5;
    ops.numThreads = 8;
    ops.firstIterMethod = 'SPA';
    ops.iterStartAt = 'heuristic';
    ops.dualGapThreshold = 0.5;

    % L = size(X, 2);
    % [spaLambdaHat, ~] = runBySPA(X, N);
    % initC = sparse(L, L);
    % [CSpaHat, ~] = computeApproxError(X, spaLambdaHat, 'maxIters', 100);
    % initC(spaLambdaHat, :) = CSpaHat;
    % ops.initC = initC;

    [hatC, Tracking] = fw(X, ops);
    [~, lambdaHat] = maxk(vecnorm(hatC, Inf, 2), N, 1);
    Tracking.hatC = hatC;
    Tracking.lambdaHatOld = lambdaHat;

end
