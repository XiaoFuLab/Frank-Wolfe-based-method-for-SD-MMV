%% Clear all things
clc; clear; close all; path(pathdef);

addpath('./fw_core/')
addpath('./utils/')

M=50; N=500; K=10; SNR=25;

%% Generate date
W = rand(M, K);
H = zeros(K, N);
H(:, 1:K) = eye(K);
H(:, K+1:end) = dirichlet_rnd(ones(K, 1), N-K);
Y = W*H;
SNR = 10^(SNR/10);
noise = randn(size(Y)); 
sigma2 = sum(vecnorm(Y, 2, 1).^2) / M / N / SNR;
noise = sqrt(sigma2)*noise;
X = Y + noise;
indices = randperm(N);
X = X(:, indices);
H = H(:, indices);
set_K_shuffled = [];
set_K = 1:K;
for ii=1:numel(set_K)
    set_K_shuffled(end+1) = find(indices == set_K(ii));
end
set_K = set_K_shuffled;

tic
%% Setting for MERIT
options = struct;
options.maxIters = 50;
options.verbose = true;
options.backend = 'matlab';
options.epsilon = 1e-6;
options.lambda = 0.5;
[C_hat_matlab, Tracking] = MERIT(X, options);
matlab_runtime = toc;

tic
options.backend = 'mex';
options.numThreads = 4;
[C_hat_mex, Tracking] = MERIT(X, options);
mex_runtime = toc;

fprintf('--------------------\n\n\n')
fprintf('Runtime of matlab version: %.2f (s) \n', matlab_runtime);
fprintf('Runtime of mex version: %.2f (s) \n', mex_runtime);
fprintf('Difference in solution between 2 versions (expected to be very small): %e\n', norm(C_hat_matlab - C_hat_mex, 'fro'));
