%% Clear all things
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')


addpath('./../')
addpath('./../utils/PrecondSPA')
addpath('./../utils')
addpath('./baseline/GeoNMF')
addpath('./baseline/SPOC')
addpath('./utils')
addpath('./../baselines/FGNSR/matlab')
addpath('./baseline/download/')

% From command line
% methodName = 'FW';
% dataSource = 'DBLP';
% dataIndex = 2;
% isCalibration = 0;
% codeName = 'something';
isSecondOrder = 0;


ExpSetting = struct;
ExpSetting.outputDir = sprintf('./results/%s/%s/%s/%d/', codeName, methodName, ...
    dataSource, dataIndex);

DataSetting = struct;
DataSetting.source = dataSource;
DataSetting.index = dataIndex;
DataSetting.isSecondOrder = isSecondOrder;

MethodSetting = struct;
MethodSetting.usingG = false;

if strcmp(methodName, 'FW') || strcmp(methodName, 'FW_200')
    MethodSetting.func = @(X, N) runByFw__tunning(X, N, lambda);
else
    error('Wrong method name')
end

MethodSetting.name = methodName;
Tracking = 0;
if ~isCalibration
    Tracking = expSetup(DataSetting, ExpSetting, MethodSetting);
end

