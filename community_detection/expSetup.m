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
    if ~isfield(DataSetting, 'source'), error('Missing field'); end;
    if ~isfield(DataSetting, 'index'), error('Missing field'); end;

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

    fprintf('Ready to run with these settings: \n');
    display(DataSetting)
    display(ExpSetting)
    display(MethodSetting)
    fprintf('----------------------------------------\n')
    %% PREPARE DATA
    [X, GT, G, Addition] = prepareData(DataSetting.source, DataSetting.index, DataSetting);
    L = size(X, 1);
    N = size(GT, 2);
    Tracking.G = G;

    fprintf('#Nodes: %d - #communities: %d\n', L, N);

    innerTic = tic;
    % if MethodSetting.usingG
    %     [thetaHat, TrackingMethod] = MethodSetting.func(X, N, G);
    % else
    [thetaHat, TrackingMethod] = MethodSetting.func(X, N);
    % end
    Tracking.thetaHat = thetaHat;
    Tracking.elapsedTime = toc(innerTic);
    Tracking.trackingMethod = TrackingMethod;
    fprintf('Elapsed time:  %.2f s\n', Tracking.elapsedTime);

    fprintf('Saved result at %s/result.mat\n', ExpSetting.outputDir);
    save(sprintf('%s/result.mat', ExpSetting.outputDir), 'Tracking');
end
