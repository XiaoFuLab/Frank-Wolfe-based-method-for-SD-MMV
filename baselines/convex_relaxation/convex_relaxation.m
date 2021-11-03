function [lambda, C] = convex_relaxation(X, epsilon, varargin) 
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
    p.addOptional('threshold', 0.05);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    L = size(X, 2);
    lambda = 2*epsilon;
    cvx_begin
        q = 2;
        variable C(L, L)
        % mix_norm = 0;
        %for i=1:size(C, 1)
        %    mix_norm = mix_norm + norm(C(i, :), q);
        %end
        minimize(sum(norms(C, q, 2)))
        % minimize(mix_norm)
        subject to
            norms(X - X*C, 2, 1) <= lambda * ones(1, L);
            zeros(size(C)) <= C;
            ones(1, L)*C == ones(1, L);
    cvx_end

    lambda = find(vecnorm(C, 2, 2)>options.threshold);
end
