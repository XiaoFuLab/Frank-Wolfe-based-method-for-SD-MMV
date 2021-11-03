function [hatA, Tracking] = estimateA2(X, S, varargin) 
% -----------------------------------------------------------------------
% estimateS Solving min_A || X - AS ||_F^2
%         st. 0 <= A 
% This is equivalent to 
%           min_A^T ||S^T A^T - X^T||_F^2
%           st. 0 <= A
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('maxIters', 10);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    % Solve || X' - S'*A'||_F^2
    [M, L] = size(X);
    N = size(S, 1);
    hatAt = rand(N, M);
    hatAt = nnlsHALSupdt(X', S', hatAt, options.maxIters);
    hatA = hatAt';
    Tracking = struct;
    Tracking.relErr = norm(X' - S'*hatAt, 'fro') / norm(X', 'fro');
end

