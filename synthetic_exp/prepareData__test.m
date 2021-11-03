%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')

M = 20; N = 10; L = 30;
[X1, pureSet, S]  = prepareData('UNIFORM', 'DIRICHLET', M, N, L, 'SNR', 1000);
S(:, pureSet)
assert(numel(pureSet) == N)

