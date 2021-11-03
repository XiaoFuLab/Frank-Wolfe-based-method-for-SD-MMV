function [hatA, hatS] = estimateCtype4(Data, K, varargin) 
% ----------------------------------------------------------------------
% 
% Given a unnormalized X, this is quite tricky way of finding A, S.
% I dont understand it

% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    X = Data.unnormalizedX;
    hatA = X(:, K);

    [hatS, Tr] = estimateS(X, hatA, ...
        'verbose', true, ...
        'maxIters', 100); % This is C' in P formulation
end
