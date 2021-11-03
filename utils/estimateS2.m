function [hatS, Tracking] = estimateS2(X, hatA, varargin) 
% -----------------------------------------------------------------------
% estimateS Solving min_S || X - AS ||_F^2
%         st. 0 <= S <= 1, 1^T S = 1^T
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('maxIters', 10); %todo check it
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    
    % Find S with NNLS 
    initS = (hatA'*hatA) \ (hatA' * X);
    hatS = nnlsHALSupdt(X, hatA, initS, 500);

    % Normalize so that S has row stochastic, as its physical interpretation
    hatS = hatS ./ (sum(hatS, 1)+eps);
end

