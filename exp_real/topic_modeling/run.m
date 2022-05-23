%% Clear all things
clc; clear; close all; path(pathdef);
addpath('../../baselines/PrecondSPA/')
addpath('../../baselines/FGNSR/matlab/')
addpath('../../utils')
addpath('../../fw_core')
addpath('../../')
addpath('./metrics/')
addpath('./utils/')

dataSource = 'tdt2';
numTopicsInput = 4;
trials = [1:10]; % Can be set to 50


for trialInd=trials
    [X, groundTruth, unnormalizedX, originalData] = prepareData(dataSource, numTopicsInput, trialInd);
    X = full(X);
    unnormalizedX = full(unnormalizedX);
    [numDocs, numVocabs] = size(X);

    fprintf(sprintf('Trial Index: %d\n', trialInd))
    fprintf(sprintf('#docs: %d - #vocab: %d\n', numDocs, numVocabs));

    [K_hat, C_hat, Tracking] = runByMERIT(X, numTopicsInput);

    Data = struct;
    Data.unnormalizedX = unnormalizedX;
    Data.originalData = originalData;
    [W_hat, H_hat] = estimate_W_H(Data, numTopicsInput, 'Addition', Tracking);

    kmeans_acc = [];
    for kmeans_trial=1:10
        pred = kmeans(W_hat, numTopicsInput, 'emptyaction', 'singleton', 'replicate', 10);
        sp = bestMap(groundTruth, pred);
        kmeans_acc(kmeans_trial) = sum(groundTruth==sp)/length(groundTruth);
    end
    acc(trialInd) = mean(kmeans_acc);

    coh(trialInd) = Coherence(H_hat', 20, originalData);
    simCount(trialInd) = sum(sum(UniequeWordsCount(H_hat', 20)))/2;
end

mean(acc)
mean(coh)
mean(simCount)

