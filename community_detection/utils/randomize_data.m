%% Clear all things
clc; clear; close all; path(pathdef);

% [X, GT] = prepareData('EMAIL', 0);
% numTrials = 10;
% for numCommunities=[10 20 30 40]
%     for trialInd=1:numTrials
%         numCommunities
%         [Xr, GTr] = random_select(X, GT, numCommunities);
%         outputDir = sprintf('%s/%s/%d/',  ...
%             './../data/email-eu-core/random-split', ...
%             num2str(numCommunities), trialInd);
%         if ~isdir(outputDir)
%             mkdir(outputDir)
%             fprintf('Created directory at: %s\n', outputDir);
%         end
%
%         size(Xr)
%         file2Save = sprintf('%s/data.mat', outputDir);
%         save(file2Save, 'Xr', 'GTr');
%     end
% end

% X = [];
% GT = [];
% for i=1:5
%     [Xs, GTs] = prepareData('DBLP', i);
%     X = blkdiag(X, Xs);
%     GT = blkdiag(GT, GTs);
% end
%
% candidates = [5 10 15 19];
% numTrials = 5;
% for i=1:numel(candidates)
%     for trialInd=1:numTrials
%         numCommunities = candidates(i);
%         [Xr, GTr] = random_select(X, GT, numCommunities);
%         outputDir = sprintf('%s/%s/%d/',  ...
%             './../data/coauthorship/random-split', ...
%             num2str(numCommunities), trialInd);
%         if ~isdir(outputDir)
%             mkdir(outputDir)
%             fprintf('Created directory at: %s\n', outputDir);
%         end
%
%         size(Xr)
%         file2Save = sprintf('%s/data.mat', outputDir);
%         save(file2Save, 'Xr', 'GTr');
%     end
% end

Data = load('./../data/youtube/induced_graph.mat');
X = Data.A;
GT = Data.GT;

candidates = [20 30 40 50 60 70 80 90 100];
numTrials = 10;
for i=1:numel(candidates)
    for trialInd=1:numTrials
        numCommunities = candidates(i);
        [Xr, GTr] = random_select(X, GT, numCommunities);
        outputDir = sprintf('%s/%s/%d/',  ...
            './../data/youtube/random-split', ...
            num2str(numCommunities), trialInd);
        if ~isdir(outputDir)
            mkdir(outputDir)
            fprintf('Created directory at: %s\n', outputDir);
        end

        size(Xr)
        file2Save = sprintf('%s/data.mat', outputDir);
        save(file2Save, 'Xr', 'GTr');
    end
end


function [X, GT] = random_select(X, GT, numCommunities, varargin) 
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
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    N = size(GT, 2);
    L = size(X, 1);
    communities = randsample(N, numCommunities);
    ind = [];
    for i=1:L
        if sum(GT(i, communities)) > 0
            ind(end+1) = i;
        end
    end
    X = X(ind, ind);
    GT = GT(ind, communities);
    GT = GT ./ sum(GT, 2);
    assert(sum(abs(sum(GT, 2) - 1) > 0.01) == 0)
end
