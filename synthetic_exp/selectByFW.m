function [lambdaHat, Tracking] = selectByFW(X, N, varargin) 
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
    p.addOptional('numOutliers', 0)
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    % FRANK-WOLFE
    ops = struct;
    ops.maxIters = 50;
    ops.verbose = 1;
    ops.debug = 0;
    ops.lambda = 0;
    ops.N = N;
    ops.tol = 0.000;
    ops.objTol = 0;
    ops.backend = 'mex';
    ops.epsilon = 1e-8;
    ops.innerVerboseInterval = 5000;
    ops.iterStartAt = 1;
    ops.firstIterMethod = 'NONE';
    
    [hatC, Tracking] = fw(X, ops);
    [v, lambdaHat] = maxk(vecnorm(hatC, Inf, 2), N+options.numOutliers);

    Tracking.hatC = hatC;
end

