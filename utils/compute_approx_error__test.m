%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common/FW')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')

d = 30; L = 500; k = 10;
d = 5; L = 20; k = 3;
A = rand(d, k);
S = drchrnd(ones(1, k), L);

pure_pixel_set = randsample(L, k);
pure_pixel_set = [1:k];
for i=1:k
    S(:, pure_pixel_set(i)) = zeros(k, 1);
    S(i, pure_pixel_set(i)) = 1;
end

X = A*S;
SNR = 10; SNR = 10^(SNR/10);
noise = randn(size(X)); 
sigma2 = sum(vecnorm(A*S, 2, 1).^2) / d / L / SNR;
noise = sqrt(sigma2)*noise;
% X = X +noise;

a = compute_approx_error(X, pure_pixel_set)

