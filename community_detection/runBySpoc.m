function [thetaHat, Tracking] = runBySpoc(X, N, isPruned, varargin) 
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
    p.addOptional('maxIters', 50);
    p.addOptional('lambda', -2);
    p.addOptional('isSparse', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    ops = p.Results;
    if isPruned
        [thetaHat, B, T] = SPOC_pruned(X, N);
    else
        [thetaHat, B, T] = SPOC(X, N);
    end
    thetaHat = thetaHat;
    Tracking = struct;
    Tracking.trackingMethod = T;
end
