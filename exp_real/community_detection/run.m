%% Clear all things
clc; clear; close all; path(pathdef);

addpath('../../utils')
addpath('../../fw_core')
addpath('../../')
addpath('utils/')

[X, GT] = prepareData('MAG', 2);
theta_hat = runByMERIT(X, size(GT, 2));

src = computeSRC(theta_hat, GT);
fprintf('SRC: %f \n', src);
