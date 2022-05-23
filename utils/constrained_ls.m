function [X] = constrained_ls(A, B, varargin) 
% ----------------------------------------------------------------------
% Use FISTA to solve the following problem
% minimize_{X}  ||AX - B||_F^2
% subject to    X \geq 0, 1^T X = 1^T

% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    % Check the options structure.
    p = inputParser;
    p.addOptional('max_iters', 100);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    ops = struct;
    ops.verbose = false;
    ops.debug = false;
    ops.max_iters = options.max_iters;
    ops.tol = 1e-10;
    p_fn = @(X, l) proj_simplex_matrix(X);
    step_size = 1/eigs(A'*A, 1);

    init_point = rand(size(A, 2), size(B, 2));
    init_point = init_point./sum(init_point);
    g_fn = @(X) A'*(A*X - B);
    ops.f_fn = @(X) norm(A*X - B)^2;
    ops.ground_truth = ones(size(init_point));
    [X, ~] = pgd_fista(g_fn, p_fn, step_size, init_point, ops);
end

