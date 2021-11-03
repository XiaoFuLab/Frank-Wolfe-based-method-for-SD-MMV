function [thetaHat, Tracking] = runBySvmCone(A, N, varargin) 
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

    [thetaHat,~] = SVM_cone_DCMMSB(A, N);
    
    Tracking = struct;
end

