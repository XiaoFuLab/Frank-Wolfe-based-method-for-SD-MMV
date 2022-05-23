%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common/FW')
addpath('~/code/matlab/common/')
addpath('~/code/matlab/common/prox_ops')
addpath('~/code/matlab/common/PGD')

addpath('./../../')

d = 30; L = 1000; k = 10;
A = rand(d, k);
S = dirichlet_rnd(ones(1, k), L);

pure_pixel_set = randsample(L, k)';
for i=1:k
    S(:, pure_pixel_set(i)) = zeros(k, 1);
    S(i, pure_pixel_set(i)) = 1;
end

X = A*S;
SNR = 40; SNR = 10^(SNR/10);
noise = randn(size(X)); 
sigma2 = sum(vecnorm(A*S, 2, 1).^2) / d / L / SNR;
noise = sqrt(sigma2)*noise;
% X = X + noise;
[C, K, ~, ~, ~, ~, t] = fgnsr(X, k, 'maxiter', 200, 'debug', 0);

disp('True basic indices: ')
disp(sort(pure_pixel_set(:))');
disp('Computed indices:')
disp(sort(K(:))');
%
% disp('==========================');
% disp('done.');
%
% figure();
% semilogy(1:numel(t.nnz), t.nnz);
% xlabel('iter')
% ylabel('NNZ')
%
