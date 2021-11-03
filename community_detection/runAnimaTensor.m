% Clear all things
clc; clear; close all; path(pathdef);

addpath('~/code/matlab/common')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')
addpath('~/code/matlab/common/tensorlab')


addpath('./../')
addpath('./../utils/PrecondSPA')
addpath('./../utils')
addpath('./baseline/GeoNMF')
addpath('./baseline/SPOC')
addpath('./utils')
addpath('./../baselines/FGNSR/matlab')
addpath('./baseline/download/')
addpath('./baseline/download/new/')
addpath('./../fw_core')

% From command line
methodName = 'ANIMA';
dataSource = 'DBLP';
dataIndex = 1;
result = [];
for dataIndex=1:5
    codeName = 'official';
    trialInd = 0;

    ExpSetting = struct;
    ExpSetting.outputDir = sprintf('./results/%s/%s/%s/%d/%d', codeName, methodName, ...
        dataSource, dataIndex, trialInd);

    DataSetting = struct;
    DataSetting.source = dataSource;
    DataSetting.index = dataIndex;
    DataSetting.isSecondOrder = 0;
    DataSetting.trialInd = trialInd;

    MethodSetting = struct;
    MethodSetting.usingG = true;
    MethodSetting.func = @(X, N) runByAnimaTensor(X, N);
    MethodSetting.name = methodName;

    for i=1:10
        Tracking = expSetup(DataSetting, ExpSetting, MethodSetting);
        DataSetting = Tracking.DataSetting;
        MethodSetting = Tracking.MethodSetting;
        DataSetting.isSecondOrder = 0;

        [X, GT] = prepareData(DataSetting.source, DataSetting.index, DataSetting);
        GT = GT(Tracking.trackingMethod.G, :);

        result(dataIndex, end+1) = computeSRC(Tracking.thetaHat, GT);
        i
    end
end
result
