function [thetaHat, Tracking] = runByAnimaTensor(A, F, varargin) 
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

    
    G = randperm(size(A, 1), 1000);
    thetaHat = anima_tensor( A, G, F );
    thetaHat = thetaHat';

    Tracking = struct;
    Tracking.G = G;
end
