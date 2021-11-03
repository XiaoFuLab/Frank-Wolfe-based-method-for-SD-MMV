%% Clear all things
clc; clear; close all; path(pathdef);
addpath('./utils')
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')
addpath('./../utils/PrecondSPA')
addpath('./../utils')
addpath('./baseline/GeoNMF')

result = [];
for index=1:3
    path2File = sprintf('./results/something/FW2order/DBLP/%d/result.mat', index);
    Tracking = expEval(path2File);
    result(index, 1) = Tracking.out;

    % N = size(Tracking.thetaHat, 2);
    % [X, GT] = prepareData('DBLP', index, 'isSecondOrder', 1);
    %
    % lambdaHat = Tracking.trackingMethod.lambdaHat;
    % hatA = X(:, lambdaHat);
    % hatS = estimateS2(X, hatA);
    % thetaHat = hatS';
    % result(index, 2) = computeSRC(thetaHat, GT);
    % Tracking.elapsedTime

    % hatA = X(:, lambdaHat);
    % thetaHat = hatA\X;
    % thetaHat = thetaHat';
    % thetaHat = max(0, thetaHat);
    % thetaHat = thetaHat ./ (sum(thetaHat, 2)+eps);
    % result(index, 3) = computeSRC(thetaHat, GT);
end
result

% Average Spearman RC: 0.079255
%
% ans =
%
%     0.1974
%
%
% ans =
%
%     0.1663
%Average Spearman RC: 0.181725
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  4.879402e-20.
% > In estimateS2 (line 19)
% In expEval_fw__run (line 22)
%
%
% ans =
%
%     0.0271
%
%
% ans =
%
%     0.0547
%
%Average Spearman RC: 0.274112
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.307923e-24.
% > In estimateS2 (line 19)
% In expEval_fw__run (line 22)
%
%
% ans =
%
%     0.1246
%
%
% ans =
%
%     0.1181
%
%
%
% %result(:,:,1) =
%
%     0.1748    0.1892    0.1892
%     0.2741    0.2698    0.1685
%     0.2497    0.2698    0.2698
%     0.2380    0.2698    0.2698
%
