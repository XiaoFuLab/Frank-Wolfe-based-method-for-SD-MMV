function [hatA, hatS] = estimateCtype3(Data, K, varargin) 
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

    % Find S with NNLS 
    initS = (hatA'*hatA) \ (hatA' * X);
    hatS = nnlsHALSupdt(X, hatA, initS, 500);

    % Normalize so that S has row stochastic, as its physical interpretation
    hatS = hatS ./ sum(hatS, 2);

    % Use this estimated S to refine A: min_A || S'*A' - X' ||_F^2
    initAt = (hatS * hatS') \ (hatS * X');
    hatAt = nnlsHALSupdt(X', hatS', initAt);
    hatA = hatAt';
    hatA = hatA ./ (vecnorm(hatA, 2, 2)+eps);
end
