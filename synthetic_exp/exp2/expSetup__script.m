%% Clear all things
% clc; clear; close all; path(pathdef);

addpath('~/code/matlab/common')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')

addpath('./../../')
addpath('./../../utils/PrecondSPA')
addpath('./../../utils')
addpath('./../../baselines/greedy_pursuit')
addpath('./../../baselines/FGNSR/matlab')
addpath('./../../fw_core')
addpath('./../')

% CONFIGS
% M=50; N=40; L=4000; SNR=10;
% methodName = 'FG';
% numTrials = 1;
% expCodeName = 'debug';
% isCalibration = 0;
%
ExpSetting = struct;
ExpSetting.outputDir = sprintf('./results/%s/%s/', expCodeName, methodName);
ExpSetting.numTrials = numTrials;

DataSetting = struct;
DataSetting.typeS = 'DIRICHLET';
DataSetting.typeA = 'UNIFORM';
DataSetting.M = M;
DataSetting.N = N;
DataSetting.L = L;
DataSetting.SNR = SNR;
DataSetting.numOutliers = 0;

MethodSetting = struct;
if strcmp(methodName, 'FW')
    MethodSetting.func = @(X, N) selectByFW(X, N, 'verbose', 0, 'numOutliers', 0);
elseif strcmp(methodName, 'FG')
    MethodSetting.func = @(X, N) selectByFastGradient(X, N);
else
    error('Mistype method name')
end

MethodSetting.name = methodName;

if ~isCalibration
    Tracking = expSetup(DataSetting, ExpSetting, MethodSetting, 'verbose', 1);
    fprintf('Succesful rate: %f \n', mean(Tracking.isSuccessful));
end

