function [hatA, Tracking] = estimateA(X, S, varargin) 
% -----------------------------------------------------------------------
% estimateS Solving min_A || X - AS ||_F^2
%         st. 0 <= A 
% This is equivalent to 
%           min_A^T ||S^T A^T - X^T||_F^2
%           st. 0 <= A
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('maxIters', 10);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    ops = struct;
    ops.verbose = false;
    ops.debug = false;
    ops.maxIters = options.maxIters;
    ops.tol = 1e-10;
    pFn = @(A, l) projNonnegative(A);
    stepSize = 1/eigs(S*S', 1);

    [D, V] = size(X);
    N = size(S, 1);
    initPoint = randn(N, D);
    gFn = @(At) S*(S'*At - X');
    ops.fFn = @(At) norm(S'*At - X', 'fro')^2;
    ops.groundTruth = sparse(size(initPoint));
    [hatAt, Tracking] = pgdFista(gFn, pFn, stepSize, initPoint, ops);
    hatA = hatAt';
    Tracking.relErr = norm(X - hatA*S, 'fro') / norm(X, 'fro');
end

