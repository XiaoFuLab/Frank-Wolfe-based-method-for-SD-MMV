%% Clear all things
clc; clear; close all; path(pathdef);
addpath('./fw_core/')
addpath('./utils/')

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

M=50; N=100; K=20; SNR=25;

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
options.lambda = 0.5;
options.epsilon = 1e-6;

options.maxIters = 200;
options.verbose = true;
options.numThreads = 4;
options.backend = 'matlab';
[C_hat, Tracking] = MERIT(X, options);
toc

fprintf('Set K: \t   ');
display(sort(set_K))
[v, set_K_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);
fprintf('Set K_hat: ')
display(sort(set_K_hat)')

figure()
tmp = xline(set_K, 'r--', 'LineWidth', 1.0, 'DisplayName', '$n \in \mathcal{K}$', 'Interpreter',  'latex')
hold on
stem(vecnorm(C_hat, Inf, 2), 'LineWidth', 1.45, 'MarkerSize', 7)
xlabel('n')
ylabel('||C(n, :)||_{\infty}')
title('Recoverying $\mathcal{K}', 'Interpreter', 'Latex')
legend(tmp(1), 'Interpreter',  'Latex')
exportgraphics(gcf, 'figure1.png', 'resolution', 300);


figure()
plot(Tracking.InnerTracking.nnz/N/N)
xlabel('Iter')
ylabel('NNZ Ratio of \textbf{C}', 'Interpreter', 'Latex')
title('nnz(C)/prod(size(C))')
exportgraphics(gcf, 'figure2.png', 'resolution', 300);

