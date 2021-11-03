%% Clear all things
clc; clear; close all; path(pathdef);
addpath('./../../synthetic_exp/')
addpath('~/code/matlab/common/FW')
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')
% M = 255; L = 10000; N = 50;
M =5; L=10; N=3;
[X, pure_pixel_set, S] = generate_data(1, 3, M, N, L, 'SNR', 10);
C = sprand(L, L, 0.2);
% C = sparse(L, L);
save('debug', 'X', 'L', 'C');
load('debug');

tic
runMultiThread(X, C);
toc

% TEST 1
% sumCe = sum(exp(C./0.01), 2);
% fprintf('\n')
% fprintf('From MATLAB \n')
% fprintf('sumC = ')
% for i=1:L
%     fprintf('%f ', full(sumCe(i)));
% end
% fprintf('\n')
%
% TEST 2
% G1 = 2*X'*(X*C - X);
% G2 = exp(C/0.01) ./ sumCe;
% G = G1 + G2;
% G'
%
% TEST 3
% tic
% for k=1:size(C, 2)
%     % sumCe(k) = full(sum(exp(C(k, :)./epsilon)));
%     X*C(:, k);
% end
% toc

% TEST 4
%
%% Clear all things
clear; close all; path(pathdef);
load('debug');
options = struct;
options.maxIters = 1;
options.objTol = -1;
options.verbose = 0;
options.lambda =1;
tic
hatC = fwCoreMatlab(X, C, 1, 1, 0.01, options);
toc
full(hatC)'
