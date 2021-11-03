function [hatS, Tracking] = estimateS(X, A, varargin) 
% -----------------------------------------------------------------------
% estimateS Solving min_S || X - AS ||_F^2
%         st. 0 <= S <= 1, 1^T S = 1^T
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('maxIters', 50); %todo check it
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    ops = struct;
    ops.verbose = options.verbose;
    ops.debug = false;
    ops.maxIters = options.maxIters;
    ops.tol = 1e-10;
    pFn = @(S, l) projSimplexMatrix(S);
    stepSize = 1/eigs(A'*A, 1);

    [M, L] = size(X);
    N = size(A, 2);

    % initPoint = randn(N, L);
    % initPoint = (A'*A)\(A'*X);
    % if isnan(sum(initPoint, 'all'))
    %     initPoint = rand(size(A, 2), size(X, 2));
    % end
    initPoint = rand(size(A, 2), size(X, 2));

    gFn = @(S) A'*(A*S - X);
    ops.fFn = @(S) norm(X - A*S, 'fro')^2;
    ops.groundTruth = sparse(size(initPoint));
    [hatS, Tracking] = pgdFista(gFn, pFn, stepSize, initPoint, ops);
    Tracking.relErr = norm(X - A*hatS, 'fro') / norm(X, 'fro');
end

