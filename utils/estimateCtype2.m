function [hatA, hatS] = estimateCtype2(Data, K, varargin) 
% ----------------------------------------------------------------------
% 
% Given a unnormalized X, this is quite tricky way of finding A, S.
% I dont understand it
% Take it from ...
% The best for C though

% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    unnormalizedX = Data.unnormalizedX;
    W = unnormalizedX(:, K);
    Ct = (W'*W)\(W'*unnormalizedX);
    Ct = nnlsHALSupdt(unnormalizedX,W,Ct,5000);
    A = Ct';
    A = bsxfun(@rdivide,A,(sum(A))+eps);

    W = (A'*A)\(A'*unnormalizedX'); % estimate W again
    W = nnlsHALSupdt(unnormalizedX',A,W);
    W = bsxfun(@rdivide,W,sqrt(sum(W.^2))+eps);
    W = full(W);
    hatA = W';
    hatS = A';
    if options.verbose
        fprintf('Relative Fitting error in finding S: %f\n', Tr.relErr);
    end
end
