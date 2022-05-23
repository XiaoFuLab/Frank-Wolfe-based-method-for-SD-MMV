function [out] = computeSRC(thetaHat, GT, varargin) 
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

    K = size(thetaHat, 2);
    R = corr(thetaHat, GT, 'Type', 'Spearman');
    M = matchpairs(-R, 1000);
    ind = sub2ind(size(R), M(:, 1), M(:, 2));
    out = mean(R(ind));
end
