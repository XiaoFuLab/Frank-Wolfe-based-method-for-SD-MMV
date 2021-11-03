%% Clear all things
clc; clear; close all; path(pathdef);
addpath('./utils')
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')
addpath('./../utils/PrecondSPA')
addpath('./../utils')
addpath('./baseline/GeoNMF')

method = 'FW';
dataset = 'MAG';
result = [];
acc = [];
for index=1:2
    for trialInd=0
        path2File = sprintf('./results/0/%s/%s/%d/%d/result.mat', ...
            method, dataset, index, trialInd);
        Tracking = expEval(path2File);
        result(index, trialInd+1) = Tracking.out;
        acc(index, trialInd+1) = Tracking.acc;
    end
end
result
% N = size(Tracking.thetaHat, 2);
% [X, GT] = prepareData('DBLP', index);
% [V, D] = eigs(X, N,'la');
% X = V';
%
% lambdaHat = Tracking.trackingMethod.lambdaHat;
% hatA = X(:, lambdaHat);
% hatS = estimateS2(X, hatA);
% thetaHat = hatS';
% computeSRC(thetaHat, GT)
% % Tracking.elapsedTime
%
