%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')
addpath('./../synthetic_exp')

M = 50; L = 20; N = 40;
[X, pure_pixel_set, S] = generate_data(1, 3, M, N, L, 'SNR', 15);
lambda = 1;


f = @(C) norm(X - X*C, 'fro')^2 + lambda*phiC(C, 1e-3);
s = @(C) exp(C/1e-3)./sum(exp(C/1e-3), 2);
g = @(C) 2*X'*(X*C - X) + lambda*s(C);
linear_app = @(C1, C2) f(C1) + trace(g(C1)'*(C2-C1));


C1 = dirichletRand(ones(L, 1), L);

C2 = dirichletRand(ones(L, 1), L);

alpha = 0.1;
C4 = C1 + alpha*(C2 -C1);
C3 = C1;
f(C4) - linear_app(C3, C4) - 0.5*5*alpha^2

