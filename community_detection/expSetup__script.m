%% Clear all things
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')


addpath('./../')
addpath('./../fw_core')
addpath('./../utils/PrecondSPA')
addpath('./../utils')
addpath('./baseline/GeoNMF')
addpath('./baseline/SPOC')
addpath('./utils')
addpath('./../baselines/FGNSR/matlab')
addpath('./baseline/download/')
addpath('./baseline/download/new/')

% From command line
% methodName = 'FW';
% dataSource = 'DBLP';
% dataIndex = 2;
% trialInd = 0;
% isCalibration = 0;
% codeName = 'something';
% isSecondOrder = 0;
% energyReserve = 0.9;


ExpSetting = struct;
ExpSetting.outputDir = sprintf('./results/%s/%s/%s/%d/%d', codeName, ...
    methodName, dataSource, dataIndex, trialInd);

DataSetting = struct;
DataSetting.source = dataSource;
DataSetting.index = dataIndex;
DataSetting.isSecondOrder = isSecondOrder;
DataSetting.trialInd = trialInd;
DataSetting.energyReserve = energyReserve;

MethodSetting = struct;
MethodSetting.usingG = false;
if strcmp(methodName, 'GEONMF')
    MethodSetting.func = @(X, N) geoNmfMethod(X, N);
elseif strcmp(methodName, 'SPOC')
    MethodSetting.func = @(X, N) runBySpoc(X, N, false);
elseif strcmp(methodName, 'SPOC_pruned')
    MethodSetting.func = @(X, N) runBySpoc(X, N, true);
elseif strcmp(methodName, 'FG')
    MethodSetting.func = @(X, N) runByFg(X, N, false);
elseif strcmp(methodName, 'FG_pruned')
    MethodSetting.func = @(X, N) runByFg(X, N, true);
elseif strcmp(methodName, 'FW') 
    MethodSetting.func = @(X, N) runByFw(X, N, false);
% elseif strcmp(methodName, 'SPACL1')
%     addpath('./baseline/download/new/')
%     MethodSetting.func = @(X, N) runBySpacl1(X, N);
elseif strcmp(methodName, 'SPACL2')
    addpath('./baseline/download/new/')
    MethodSetting.func = @(X, N) runBySpacl2(X, N, false);
elseif strcmp(methodName, 'SPACL2_pruned')
    addpath('./baseline/download/new/')
    MethodSetting.func = @(X, N) runBySpacl2(X, N, true);
% elseif strcmp(methodName, 'SVMCONE')
%     addpath('./baseline/download/new/SVM_cone')
%     MethodSetting.func = @(X, N) runBySvmCone(X, N);
else
    error('Wrong method name')
end

MethodSetting.name = methodName;
Tracking = 0;
if ~isCalibration
    Tracking = expSetup(DataSetting, ExpSetting, MethodSetting);
end

