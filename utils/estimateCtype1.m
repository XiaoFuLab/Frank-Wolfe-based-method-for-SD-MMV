function [hatA, hatS] = estimateCtype1(Data, K, varargin) 
% ----------------------------------------------------------------------
% 
% Given a normalized X, this is the most straighforward way to find A, S.


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

    [hatS, Tr] = estimateS(X, hatA, 'verbose', true, 'maxIters', 100); 
    
    if options.verbose
        fprintf('Relative Fitting error in finding S: %f\n', Tr.relErr);
    end
end
