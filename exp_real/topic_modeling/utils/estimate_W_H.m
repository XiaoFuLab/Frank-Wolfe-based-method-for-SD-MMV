function [W_hat, H_hat] = estimate_W_H(Data, K, varargin) 
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    originalData = Data.originalData;
    X = Data.unnormalizedX;
    W_hat = X(:, K);
    disp('sum(W_hat) is')
    disp(sum(W_hat))

    % Find S with NNLS 
    initS = (W_hat'*W_hat) \ (W_hat' * X);
    if isnan(sum(initS, 'all'))
        initS = rand(numel(K), size(X, 2));
    end

    H_hat = nnlsHALSupdt(X, W_hat, initS, 500);

    % Normalize so that S has row stochastic, 
    % as its physical interpretation
    H_hat = H_hat ./ (sum(H_hat, 2));

    % Use this estimated S to refine A: min_A || S'*A' - X' ||_F^2
    initAt = (H_hat * H_hat') \ (H_hat * originalData');
    hatAt = nnlsHALSupdt(originalData', H_hat', initAt);
    W_hat = hatAt';
    W_hat = W_hat ./ (vecnorm(W_hat, 2, 2)+eps);
end
