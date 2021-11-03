function [hatA, hatS] = estimatePtype1(Data, K, varargin) 
% ----------------------------------------------------------------------
% 


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    originalData = Data.originalData;
    X = Data.unnormalizedX;
    hatA = X(:, K);
    disp('sum(hatA) is')
    disp(sum(hatA))

    % Find S with NNLS 
    initS = (hatA'*hatA) \ (hatA' * X);
    if isnan(sum(initS, 'all'))
        initS = rand(numel(K), size(X, 2));
    end

    hatS = nnlsHALSupdt(X, hatA, initS, 500);

    % Normalize so that S has row stochastic, 
    % as its physical interpretation
    hatS = hatS ./ (sum(hatS, 2));

    % Use this estimated S to refine A: min_A || S'*A' - X' ||_F^2
    initAt = (hatS * hatS') \ (hatS * originalData');
    hatAt = nnlsHALSupdt(originalData', hatS', initAt);
    hatA = hatAt';
    hatA = hatA ./ (vecnorm(hatA, 2, 2)+eps);

    if options.verbose
        fprintf('Relative Fitting error in finding S: %f\n', Tr.relErr);
    end
end
