%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')

addpath('./baselines/PrecondSPA/')
addpath('./../')
addpath('./../baselines/greedy_pursuit')
addpath('./../baselines/FGNSR/matlab/')
addpath('./metrics/')
addpath('./utils/')
addpath('./../utils/')


Acc = [];
Coh = [];
SimCount = [];
for trialInd=49:50
    trialInd
    % outputDir = sprintf('./results/tdt2/P_PCA2/tfidf/FW_PCA_test/numTopics_5/trial_%d', trialInd);
    outputDir = sprintf('./results/tdt2/P_PCA2/tfidf/FG_official/numTopics_7/trial_%d', trialInd);
    estimateASFunc = @estimatePtype1;
    Tracking = loadIt(outputDir, estimateASFunc, 'verbose', true, 'trialIndex', trialInd);
    Acc(end+1) = Tracking.acc(trialInd);
    Coh(end+1) = Tracking.coh(trialInd);
    SimCount(end+1) = Tracking.simCount(trialInd);
    % fprintf('Acc: %f - Coh: %f - SimCount: %f \n', mean(Tracking.acc), ...
    %     mean(Tracking.coh), mean(Tracking.simCount));
end
mean(Acc)
mean(Coh)
mean(SimCount)

function [hatA, hatS] = estimatePTypeForFastAnchor(Data, K, varargin) 
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
    p.addOptional('Addition', struct);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    originalData = Data.originalData;
    X = Data.unnormalizedX;
    hatA = options.Addition.hatA;

    % Find S with NNLS 
    initS = (hatA'*hatA) \ (hatA' * X);
    if isnan(sum(initS, 'all'))
        initS = rand(numel(K), size(X, 2));
    end

    hatS = nnlsHALSupdt(X, hatA, initS, 500);

    % Normalize so that S has row stochastic, 
    % as its physical interpretation
    hatS = hatS ./ (sum(hatS, 2));

    % Use this estimated S to refine A: min_A || S'*A' - X' ||_F^2
    initAt = (hatS * hatS') \ (hatS * originalData');
    hatAt = nnlsHALSupdt(originalData', hatS', initAt);
    hatA = hatAt';
    hatA = hatA ./ (vecnorm(hatA, 2, 2)+eps);
end
