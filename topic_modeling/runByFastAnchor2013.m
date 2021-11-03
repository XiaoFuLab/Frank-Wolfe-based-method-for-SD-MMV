function [lambdaHat, hatC, Tracking] = runByFastAnchor2013(X, N, varargin) 
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

    assert(size(X, 1) == size(X, 2))
    [A, R, J] = arora(X, N, 2013);
    hatC = nan;
    lambdaHat = nan;
    Tracking = struct;
    Tracking.hatA = A;
end

