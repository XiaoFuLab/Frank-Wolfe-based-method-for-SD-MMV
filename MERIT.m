function [C, Tracking] = MERIT(X, varargin) 
% -------------------------------------------------------------------------
%% MERIT: a FW-based algorithm for solving 
%           minimize_{C}   ||X-XC||_F^2 + \lambda \varPhi_{\mu}(C)
%           subject to     C \geq 0, 1^T C = 1^T
%
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)
%
% Papams:
%   Required params
% - X: data input
%
%   Hyperparameters 
% - lambda: lambda in the objective function above
% - epsilon: hyperparameter of the regularization term (default=1e-5)
%
%   Settings for the main procedure
% - objTol: a threshold in relative change of objective function 
% - maxIters: max number of iterations
% - backend: either `matlab` or `mex`
% - numThreads: only be applicable if backend=`mex`. It enables matrix multiplication to be performed parallel.
% - estimatedN: only be applicable if backend=`mex`. It is used to allocate memory at the beginning.
% - initC: initilization for C
% - iterStartAt: standard Frank-Wolfe start with `iterStartAt`=1, although a warm version can use `iterStartAt`>0 (default 1)
% - dualGapThreshold: another criterion for early stopping (default 1e-5)
%
%   Logging behaviours
% - verbose: logging or not
% - debug: add some more information for debugging, might require more computation
% - innerVerboseInterval: just logging interval
%
% Output:
% - C: optimal C
% - Tracking: a structure with some information collecting during the process.
%
% Reference: 
% Nguyen, T., Fu, X. and Wu, R., 2021. Memory-efficient convex optimization for 
% self-dictionary separable nonnegative matrix factorization: A frank-wolfe approach. 
% arXiv preprint arXiv:2109.11135. (Submitted to TSP, accepted, 2022.)
% ------------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('lambda', 1);
    p.addOptional('epsilon', 1e-5);

    p.addOptional('objTol', 0);
    p.addOptional('maxIters', 500);
    p.addOptional('backend', '');
    p.addOptional('numThreads', 4);
    p.addOptional('estimatedN', 30);
    p.addOptional('initC', nan); 
    p.addOptional('iterStartAt', 1); 
    p.addOptional('dualGapThreshold', 1e-5); 

    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('innerVerboseInterval', 1000);

    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    [~, N] = size(X);

    Tracking = struct;
    Tracking.obj = zeros(options.maxIters, 1);
    Tracking.time = zeros(options.maxIters, 1);
    Tracking.stepSize = zeros(options.maxIters, 1);
    Tracking.endingCode = 0;
    Tracking.checkpoint = [];

    epsilon = options.epsilon;

    nnzX = nnz(X)/prod(size(X));
    if options.verbose
        fprintf('X nnz ratio: %f \n', nnzX);
    end

    XT = X';
    lambda = options.lambda;
    if lambda<0
        error('lambda must be nonnegative')
    end
    Tracking.lambda = options.lambda;

    if ~isnan(options.initC)
        C = options.initC;
    else
        C = sparse(N, N);
    end
    
    iterStart = options.iterStartAt;
    maxIters = options.maxIters;
    if maxIters>iterStart
        if strcmp(options.backend, 'matlab')
            [C, InnerTracking] = fwCoreMatlab(X, C, lambda, iterStart, options.epsilon, ...
                maxIters, options.innerVerboseInterval, ...
                options.objTol, options.verbose, options.dualGapThreshold, 'debug', options.debug);
            Tracking.InnerTracking = InnerTracking;
        elseif strcmp(options.backend, 'mex')
            estimatedN = int64(options.estimatedN);
            numThreads = options.numThreads;
            C = fwCoreMex(X, C, lambda, iterStart, options.epsilon, ...
                maxIters, options.innerVerboseInterval, ...
                options.objTol, options.verbose, ...
                estimatedN, numThreads, options.dualGapThreshold);
        else
            error('backend is either mex or matlab');
        end
    end
end

