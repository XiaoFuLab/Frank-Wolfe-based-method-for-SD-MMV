%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')

addpath('./synthetic_exp')
addpath('./topic_modeling/baselines/PrecondSPA')
addpath('./utils')

addpath('./fw_core/')

M = 380; L = 1000; N = 4;
% M = 5; L = 5; N = 3;

[X, pure_pixel_set, S] = generate_data(1, 3, M, N, L, 'SNR', 30);
save('debug.mat', 'X', 'pure_pixel_set');
load('debug.mat', 'X', 'pure_pixel_set')

options = struct;
options.maxIters = 30;
options.lambda = 10;
options.objTol = -1;
options.N = N;
options.checkpointIterval = 1;
options.firstIterMethod = 'NONE';
options.verbose = 1;
options.innerVerboseInterval = 6000;
options.numThreads = 4;

% fprintf('From MATLAB TUNNING \n');
% tic
% options.backend = 'matlabTunning';
% options.epsilon = 1e-6;
% [C_hat4, Tracking] = fw(X, options);
% toc
% fprintf('\n');
%
% fprintf('From MATLAB\n');
% tic
% options.backend = 'matlab';
% options.epsilon = 1e-6;
% [C_hat1, Tracking] = fw(X, options);
% toc
% sum(vecnorm(C_hat1, Inf, 2))
% fprintf('\n');

fprintf('From MEX dense\n');
tic
options.backend = 'mex';
options.epsilon = 1e-6;
[C_hat2, Tracking] = fw(X, options);
toc

% fprintf('From MEX sparse \n');
% tic
% options.backend = 'mex';
% options.epsilon = 1e-6;
% [C_hat3, Tracking] = fw(sparse(X), options);
% toc
%
%
% fprintf('Difference between matlab and dense mex: %e \n', norm(C_hat2 - C_hat1, 'fro'));
% fprintf('Difference between dense and sparse mex: %e \n', norm(C_hat2 - C_hat3, 'fro'));
% fprintf('Difference between matlab and sparse mex: %e \n', norm(C_hat1 - C_hat3, 'fro'));
% fprintf('Difference between matlab and matlab tunning: %e \n', norm(C_hat1 - C_hat4, 'fro'));


% fprintf('\n');
% sort(pure_pixel_set)
% [v, lambdaHat] = maxk(vecnorm(C_hat2, Inf, 2), N, 1);
% sort(lambdaHat)'

% v = vecnorm(C_hat2, 2, 2);
% figure();
% scatter(1:L, v);
% title('L2')


% v = vecnorm(C_hat1, Inf, 2);
% figure();
% scatter(1:L, v);
% title('L-Inf C-hat1')
%
% v = vecnorm(C_hat2, Inf, 2);
% figure();
% scatter(1:L, v);
% title('L-Inf C-hat2')

% [v, lambdaHat] = maxk(vecnorm(C_hat2, Inf, 2), N, 1);
% full(v)'
% lambdaHat'
%
% nnz(C_hat2) /prod(size(C_hat2))
