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
for index=5:5
        path2File = sprintf('./results/official/FW_200/DBLP/%d/result.mat', ...
            index);
        Tracking = expEval(path2File);
        result(end +1) = Tracking.out;

        N = size(Tracking.thetaHat, 2);
        [X, GT] = prepareData(Tracking.DataSetting.source, index, Tracking.DataSetting);
        [V, D] = eigs(X, N,'la');
        X = V';
        
        for j=1:numel(Tracking.trackingMethod.checkpoint)

            hatC = Tracking.trackingMethod.checkpoint(j).C;
            [~, lambdaHat] = maxk(vecnorm(hatC, 2, 2), N, 1);
            hatA = X(:, lambdaHat);
            hatS = estimateS(X, hatA, 'maxIters', 100);
            thetaHat = hatS';
            result(end+1) = computeSRC(thetaHat, GT);

            % thetaHat = hatA\V';
            % thetaHat = thetaHat';
            % thetaHat = max(0, thetaHat);
            % thetaHat = thetaHat ./ (sum(thetaHat, 2)+eps);
            % result(end+1) = computeSRC(thetaHat, GT);
        end

        %
        % theta = X'/hatA';
        % tre = 1e-12;
        % theta(theta<tre) = 0;
        % [ theta ] = normalize_row_l1_( theta );
        % thetaHat = theta;
        % result(index, 4) = computeSRC(thetaHat, GT);
        
end
result

function [ theta_l1 ] = normalize_row_l1_( theta )

d_theta = sum(theta,2);
S = find(d_theta>0);
theta_l1 = theta;
theta_l1(S,:) = bsxfun(@times, theta(S,:), 1./(sum(theta(S,:), 2)));

end

