%% Clear all things
clc; clear; close all; path(pathdef);

L = 10000;
X = rand(L, L);
C = sprand(1, L, 0.2);

tic
(C*X)*X';
toc

tic
testSparse(C, X);
toc
