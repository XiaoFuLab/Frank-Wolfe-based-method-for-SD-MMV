% 
% 


% ---------------------------------------------------------
addpath('../../baselines/PrecondSPA/')
addpath('../../baselines/FGNSR/matlab/')
addpath('../../utils')
addpath('../../fw_core')
addpath('../../')

% N=500;
% method_name='MERIT';
M=50; K=40; SNR=10;
alpha = ones(K, 1);

W = rand(M, K);
H = zeros(K, N);
H(:, 1:K) = eye(K);
H(:, K+1:end) = dirichlet_rnd(alpha, N-K);
Y = W*H;
SNR = 10^(SNR/10);
noise = randn(size(Y)); 
sigma2 = sum(vecnorm(Y, 2, 1).^2) / M / N / SNR;
noise = sqrt(sigma2)*noise;
X = Y + noise;
indices = randperm(N);
X = X(:, indices);
H = H(:, indices);
r_pure_pixel_set = [];
pure_pixel_set = 1:K;
for ii=1:numel(pure_pixel_set)
    r_pure_pixel_set(end+1) = find(indices == pure_pixel_set(ii));
end
pure_pixel_set = r_pure_pixel_set;


if strcmp(method_name, 'FastGradient')
    fg_tic = tic;
    [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
        'debug', 0, 'verbose', true);
    fprintf('Run FastGradient.... Done\n')
elseif strcmp(method_name, 'MERIT')
    Lambda_hat_SPA = SPAselect(X, K);
    fw_tic = tic;
    options = struct;
    options.epsilon = 1e-5;
    W_hat_spa = X(:, Lambda_hat_SPA);
    H_hat_spa = constrained_ls(W_hat_spa, X, 'max_iters', 100);
    %   Use SPA solution to determine lambda
    options.lambda = norm(X - W_hat_spa*H_hat_spa, 'fro')/K;

    options.maxIters = 100;
    options.verbose = 0;
    options.debug = 0;
    options.backend = 'mex';
    [C_hat, fw_tracking] = MERIT(X, options);
    [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);
    fprintf('Run MERIT.... Done\n')
else
    error('method_name must be either MERIT or FastGradient')
end
