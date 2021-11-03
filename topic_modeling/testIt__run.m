%% Clear all things
% clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')

addpath('./baselines/PrecondSPA/')
addpath('./baselines/')
addpath('./baselines/snpa_v3/')
addpath('./../')
addpath('./../baselines/greedy_pursuit')
addpath('./../baselines/FGNSR/matlab/')
addpath('./../baselines/FGNSR/matlab/')
addpath('./metrics/')
addpath('./utils/')

% methodName = 'FastAnchor';
% dataSource = 'tdt2';
% dataType = 'tfidf';
% modelType = 'P_PCA2';
% numTopicsInput = 6;
% nTrials = 1;
% trialIndex = 1;
% isCalibration = 0;

ExpSetting = struct;
ExpSetting.outputDir = sprintf('./results/%s/%s/%s/%s/numTopics_%d/', ...
    dataSource, modelType, dataType, methodName, numTopicsInput);
ExpSetting.logInterval = 1;

DataSetting = struct;
DataSetting.nTopicsArr = [numTopicsInput];
DataSetting.nTrials = nTrials;
DataSetting.source = dataSource;
DataSetting.dataType = dataType;
DataSetting.modelType = modelType;

MethodSetting = struct;

isNormalized = true;
if strcmp(methodName, 'FastAnchor2012')
    MethodSetting.func = @runByFastAnchor2012;
    DataSetting.isDense = true;
    isNormalized = false;
elseif strcmp(methodName, 'FastAnchor2013')
    MethodSetting.func = @runByFastAnchor2013;
    DataSetting.isDense = true;
    isNormalized = false;
elseif strcmp(methodName, 'SPA') || strcmp(methodName, 'SPA_PCA')
    MethodSetting.func = @runBySPA;
elseif strcmp(methodName, 'XRAY')
    DataSetting.isDense = true;
    MethodSetting.func = @runByXray;
elseif strcmp(methodName, 'LDA')
    MethodSetting.func = @runByLDA;
else
    error('Typo');
end

MethodSetting.name = methodName;
Tracking = testIt(DataSetting, ExpSetting, MethodSetting, ...
    'trialIndex', trialIndex, ...
    'isNormalized', isNormalized);

