function [acc] = computeAcc(thetaHat, GT, varargin) 
% ----------------------------------------------------------------------
% 
% This is average Spearman ranking correlation


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;


    % [~, thetaHat] = max(thetaHat, [], 2);
    % [~, GT] = max(GT, [], 2);
    % acc = mean(mean(thetaHat==GT));
    

    numCommunities = size(thetaHat, 2);
    thetaHat = thetaHat >= 1/numCommunities;
    GT = GT >= 1/numCommunities;
    acc = mean(mean(thetaHat==GT));
end
