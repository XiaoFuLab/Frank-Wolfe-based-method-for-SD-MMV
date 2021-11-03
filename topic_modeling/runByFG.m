function [lambdaHat, hatC, Tracking] = runByFG(X, N, varargin) 
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
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    [X, K, ~, ~, ~, ~, t] = fgnsr(X, N, 'maxiter', 50, ...
        'debug', 0, 'verbose', true, 'log_interval', 1, 'init', 2);
    Tracking = t;
    lambdaHat = K;
    hatC = X;
end

