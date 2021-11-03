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

M = 3; L = 5; N = 2;

[X, pure_pixel_set, S] = generate_data(1, 3, M, N, L, 'SNR', 30);
save('debug.mat', 'X', 'pure_pixel_set', 'S');
load('debug.mat', 'X', 'pure_pixel_set', 'S');
fprintf('Optimal fitting error: %f \n', norm(X - X(:, pure_pixel_set)*S, 'fro'));
X = [X X(:, pure_pixel_set(1)) ];
S = [ S S(:, pure_pixel_set(1))];
pure_pixel_set = pure_pixel_set;

options = struct;
options.maxIters = 100;
options.lambda = 0.01;
options.objTol = -1;
options.N = N;
options.checkpointIterval = 1;
options.firstIterMethod = 'NONE';
options.verbose = 1;
options.innerVerboseInterval = 6000;
options.numThreads = 4;
options.iterStartAt = 1;


fprintf('From MATLAB\n');
tic
options.backend = 'matlab';
options.epsilon = 1e-6;
[C_hat, Tracking] = fw(X, options);
toc
sum(vecnorm(C_hat, Inf, 2))
fprintf('\n');


fprintf('\n');
% disp('Pure Pixel: ')
% disp(sort(pure_pixel_set))
[v, lambdaHat] = maxk(vecnorm(C_hat, Inf, 2), N, 1);
% disp('Lambda Hat: ')
% disp(sort(lambdaHat)')
fprintf('is successfull: %d \n', all(sort(lambdaHat(:)) == sort(pure_pixel_set(:))))

v = vecnorm(C_hat, Inf, 2);
figure();
L = size(X, 2);
scatter(1:L, v, 'filled');
title(sprintf('L-Inf with lambda=%f', options.lambda))

[hatS, ~] = compute_approx_error(X, lambdaHat);
fprintf('FW fitting error: %f \n', norm(X - X(:, lambdaHat)*hatS, 'fro'));

spaLambdaHat = SPAselect(X, N);
% disp(sort(spaLambdaHat(:)'))
fprintf('is SPA successfull: %d \n', all(sort(spaLambdaHat(:))==sort(pure_pixel_set(:))))

[hatS, ~] = compute_approx_error(X, spaLambdaHat);
fprintf('FW fitting error: %f \n', norm(X - X(:, spaLambdaHat)*hatS, 'fro'));

fprintf('Optimal fitting error: %f \n', norm(X - X(:, pure_pixel_set)*S, 'fro'));
