function [out, tracking] = compute_approx_error(X, K, varargin) 
% -------------------------------------------------------------------------
% 
% Find S 
% \minimize 0.5*\norm{X - X(:, K)*H, 'fro'}^2
% Subject to each column of H is a probability simplex.
% return: minimum value / norm(X, 'fro')^2

% Require:
% - pgd_fista
% - proj_simplex
%
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -------------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('max_iters', 100);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

   
    k = numel(K);
    L = size(X, 2);
    A = X(:, K);

    ops = struct;
    ops.verbose = false;
    ops.debug = false;
    ops.max_iters = options.max_iters;
    ops.tol = 1e-10;
    ops.max_iters = options.max_iters;
    p_fn = @(H, l) proj_simplex_matrix(H);
    step_size = 1/eigs(A'*A, 1);

    init_point = randn(k, L);
    g_fn = @(H) A'*(A*H - X);
    ops.f_fn = @(H) norm(X - A*H)^2;
    ops.ground_truth = ones(size(init_point));
    [S, tracking] = pgd_fista(g_fn, p_fn, step_size, init_point, ops);
    out = S;

end
