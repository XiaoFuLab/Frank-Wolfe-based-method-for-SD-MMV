function [thetaHat, Tracking] = runBySpacl2(A, N, isPruned, varargin) 
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
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    ops = p.Results;

    if isPruned
        [thetaHat,~, t] = SPACL_pruned(A, N, 1);
    else
        [thetaHat,~, t] = SPACL(A, N, 1);
    end
    
    Tracking = struct;
    Tracking.lambdaHat = t.lambdaHat;
    Tracking.I_outlier = t.I_outlier;
end

