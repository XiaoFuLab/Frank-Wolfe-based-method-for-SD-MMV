clc; clear; close all;

A = rand(100, 10);
X = rand(10, 10);
X = X./sum(X);
B = A*X;

X_hat = constrained_ls(A, B, 'max_iters', 500);
fprintf('Distance to optimal solution: %e \n', norm (X - X_hat, 'fro'))

