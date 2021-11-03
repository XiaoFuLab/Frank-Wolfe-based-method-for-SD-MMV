function [X, I] = truncateData(X, keep, varargin) 
% ----------------------------------------------------------------------
% 
% X: L x m
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

    colEnergy = vecnorm(X, 2, 2).^2;
    sortColEnergy = sort(colEnergy, 'descend');
    totalEnergy = sum(colEnergy);
    sumEnergy = 0;
    threshold = 0;
    L = size(X, 1);
    for i=1:L
        sumEnergy = sumEnergy + sortColEnergy(i);
        if sumEnergy > keep*totalEnergy
            threshold = sortColEnergy(i);
            break
        end
    end
    zeroFiltering = find(colEnergy < threshold);
    I = zeroFiltering;

    X(zeroFiltering, :) = [];
end
