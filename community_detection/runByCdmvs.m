function [thetaHat, Tracking] = runByCdmvs(Y, N, G, varargin) 
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

    [thetaHat,~] = CDMVS(Y);
    thetaHat = thetaHat';
    
    Tracking = struct;
    Tracking.G = G;
end

