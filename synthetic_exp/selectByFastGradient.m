function [lambdaHat, t] = selectByFastGradient(X, K, varargin) 
% ----------------------------------------------------------------------
% 
% Summary of this function goes here
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

    gp_options = struct;
    gp_options.rank = K;

    [~, lambdaHat, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 50, ...
        'debug', 0, 'verbose', false);
end
