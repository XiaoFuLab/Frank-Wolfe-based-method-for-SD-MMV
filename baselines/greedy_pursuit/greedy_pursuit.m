function [lambda] = greedy_pursuit(X, varargin) 
% -------------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -------------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('rank', 5);
    p.addOptional('epsilon', 1);
    p.addOptional('p', 'Inf'); % p-Norm: 1, 2, 'Inf'
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    q = options.p;
    i = 0;
    R = X;
    delta = 2*options.epsilon;
    is_stopped = false;
    while is_stopped == false
        i = i+1;
        [~, idx] = max(vecnorm(R'*X, q, 1));
        lambda(i) = idx;

        if i>1
            if options.rank > 0
                is_stopped = i > options.rank;
            else
                is_stopped = check_criteria(X(:, idx), Lambda, delta, options);
            end
        end
        if is_stopped
            break
        end

        Lambda = X(:, lambda);
        R = ((eye(size(X, 1))) - Lambda*inv(Lambda'*Lambda)*Lambda') *X;
    end
    lambda = lambda(1:end-1);


    function [out] = check_criteria(x, X, delta, varargin) 
    % -------------------------------------------------------------------------
    %
    % Summary of this function goes here
    % Detailed explanation goes here


    % Author: Tri Nguyen (nguyetr9@oregonstate.edu)

    % -------------------------------------------------------------------------
        
        % Check the options structure.
        p = inputParser;
        p.addOptional('verbose', true);
        p.addOptional('debug', true);
        p.addOptional('max_iters', 500);
        p.addOptional('f_fn', @(c) norm(x - X*c)^2);
        p.addOptional('ground_truth', zeros(size(X, 2), 1));
        p.addOptional('prox_weight', 0);
        p.KeepUnmatched = true;
        p.parse(varargin{:});
        options = p.Results;
        
        g_fn = @(c) 2*X'*(X*c - x);
        p_fn = @proj_simplex;
        step_size = 1/eigs(X'*X, 1)/2;
        init_point = randn(size(X, 2), 1);

        [c_hat, tracking] = pgd_fista(g_fn, p_fn, step_size, init_point, ...
            options);
        if options.verbose, fprintf('Num iters: %d\n', tracking.num_iters); end;
        
        obj = norm(x - X*c_hat);
        out = obj <= delta;
    end
end
