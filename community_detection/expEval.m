function [Tracking] = expEval(path2File, varargin) 
% ----------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('special', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    dumpedData = load(path2File, 'Tracking');
    Tracking = dumpedData.Tracking;
    DataSetting = Tracking.DataSetting;
    MethodSetting = Tracking.MethodSetting;

    [X, GT] = prepareData(DataSetting.source, ...
        DataSetting.index, DataSetting);
    if MethodSetting.usingG
        GT = GT(Tracking.trackingMethod.G, :);
    end

    if ~options.special
        Tracking.out = computeSRC(Tracking.thetaHat, GT);
        Tracking.acc = computeAcc(Tracking.thetaHat, GT);
    else
        Tracking.out = 0;
    end

    display(Tracking.MethodSetting)
    display(Tracking.DataSetting)
    fprintf('Average Spearman RC: %f\n', Tracking.out);
end
