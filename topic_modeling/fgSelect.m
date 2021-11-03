function [K, X, t] = fgSelect(X, N, varargin) 
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
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    [X, K, ~, ~, ~, ~, t] = fgnsr(X, N, 'maxiter', 20, ...
        'debug', 0, 'verbose', true, 'log_interval', 1);
end
