function [Tracking] = expSetup(DataSetting, ExpSetting, MethodSetting, varargin) 
% expSetup    Go for it
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)
%
%% Inputs:
%       DataSetting : A structure contains all settings related to dataset
%           - nTopicsArr    : an 1D array contains list of number of topics
%           - nTrials       : Number of trials
%           - ...
%       ExpSetting  : ababa , and even if 
%                     it gets longer
%       MethodSetting   : A structure contains information about method to test
%           - func          : a function handler that receives 2 params: X, N,
%                           and outputs K
%           - name          : method name
%                       
%% Outputs:
%       Tracking    : A structure contains predicting result
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    if options.debug
        fprintf('WARNING: DEBUGGING IS ON!!!! \n')
        pause(2)
    end

    %% DATA SETTING
    if ~isfield(DataSetting, 'typeA'), error('Missing field'); end;
    if ~isfield(DataSetting, 'typeS'), error('Missing field'); end;
    if ~isfield(DataSetting, 'M'), error('Missing field'); end;
    if ~isfield(DataSetting, 'N'), error('Missing field'); end;
    if ~isfield(DataSetting, 'L'), error('Missing field'); end;
    if ~isfield(DataSetting, 'SNR'), error('Missing field'); end;

    %% EXPERIMENT SETTING
    if ~isfield(ExpSetting, 'outputDir'), error('Missing field'); end;

    %% METHOD SETTING
    if ~isfield(MethodSetting, 'func'), error('Missing field'); end;
    if ~isfield(MethodSetting, 'name'), error('Missing field'); end;

    Tracking = struct;
    Tracking.elapsedTime = 0;
    Tracking.DataSetting = DataSetting;
    Tracking.ExpSetting = ExpSetting;
    Tracking.MethodSetting = MethodSetting;

    if ~isdir(ExpSetting.outputDir)
        mkdir(ExpSetting.outputDir)
        fprintf('Created directory at: %s\n', ExpSetting.outputDir);
    end

    if options.verbose
        fprintf('Ready to run with these settings: \n');
        display(DataSetting)
        display(ExpSetting)
        display(MethodSetting)
        fprintf('----------------------------------------\n')
    end

    Tracking.lambdaHat = cell(ExpSetting.numTrials, 1);
    Tracking.elapsedTime = zeros(ExpSetting.numTrials, 1);
    Tracking.trackingMethod = cell(ExpSetting.numTrials, 1);
    Tracking.GT = cell(ExpSetting.numTrials, 1);
    Tracking.coverRate = zeros(ExpSetting.numTrials, 1);
    Tracking.isSuccessful = zeros(ExpSetting.numTrials, 1);

    for trialInd=1:ExpSetting.numTrials
        
        %% PREPARE DATA
        [X, GT] = prepareData(DataSetting.typeA, DataSetting.typeS, ...
            DataSetting.M, DataSetting.N, DataSetting.L, ...
            'SNR', DataSetting.SNR, DataSetting);
        L = size(X, 2);
        N = numel(GT);

        if options.verbose
            fprintf('#L: %d - #N: %d\n', L, N);
            fprintf('||X||_F = %f \n', norm(X, 'fro'));
        end

        innerTic = tic;

        [lambdaHat, TrackingMethod] = MethodSetting.func(X, N);
        Tracking.lambdaHat{trialInd} = lambdaHat;
        Tracking.elapsedTime(trialInd) = toc(innerTic);
        Tracking.trackingMethod{trialInd} = TrackingMethod;
        Tracking.GT{trialInd} = GT;
        Tracking.coverRate(trialInd) = sum(ismember(GT, lambdaHat))/...
            numel(lambdaHat);
        Tracking.isSuccessful(trialInd) = ...
            (numel(setdiff(lambdaHat, GT)) == 0) && ...
            (numel(setdiff(GT, lambdaHat)) == 0);

        if options.verbose
            fprintf('Elapsed time:  %.2f s\n', Tracking.elapsedTime);
        end
    end

    pathToFile = sprintf('%s/result.mat', ExpSetting.outputDir);
    save(pathToFile, 'Tracking');
    fprintf('Saved result at %s \n', pathToFile);
end
