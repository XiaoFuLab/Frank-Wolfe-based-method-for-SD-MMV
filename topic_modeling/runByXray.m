function [lambdaHat, hatC, Tracking] = runByXray(X, N, varargin) 
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

    [lambdaHat, ~] = SNPA(X, N);
    hatC = nan;
    Tracking = struct;
end

