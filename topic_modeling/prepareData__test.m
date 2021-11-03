%% Clear all things
clc; clear; close all; path(pathdef);
addpath('./utils/')


nTopics = 3:10;
sizes1 = [];
sizes2 = [];
lips = [];
for i=1:numel(nTopics)
    [barX, groundTruth, X] = prepareData('tdt2', nTopics(i), 1, ...
    'dataType', 'tfidf', ...
    'modelType', 'P_PCA2', ...
    'verbose', true);
    sizes1(end+1) = size(barX, 2);
    sizes2(end+1) = size(barX, 1);
    lips(end+1) = eigs(barX'*barX, 1, 'la');
end
sizes1
sizes2
return

nTopics = 2;
trialInd = 1;
[barX, groundTruth, X] = prepareData('reuters21578', nTopics, trialInd, ...
    'dataType', 'tf', ...
    'modelType', 'C');

tmp = sum(barX, 1);
assert( norm(tmp - 1, 'fro') <= 1e-10)
assert(all(size(barX) == size(X)))
[nDocs, nVocabs] = size(barX);
assert(nDocs < nVocabs)
assert(nDocs == numel(groundTruth))
assert(max(X(:)) > 1)
assert(max(barX(:)) <= 1)

[barX, groundTruth, X] = prepareData('reuters21578', nTopics, trialInd, ...
    'dataType', 'tfidf', ...
    'modelType', 'C');

tmp = sum(barX, 1);
assert( norm(tmp - 1, 'fro') <= 1e-10)
assert(all(size(barX) == size(X)))
[nDocs, nVocabs] = size(barX);
assert(nDocs < nVocabs)
assert(nDocs == numel(groundTruth))
assert(max(barX(:)) <= 1)

[barX, groundTruth, X] = prepareData('reuters21578', nTopics, trialInd, ...
    'dataType', 'tf', ...
    'modelType', 'P_CLEAN');

tmp = sum(barX, 1);
assert( norm(tmp - 1, 'fro') <= 1e-10)
assert(all(size(barX) == size(X)))
[nDocs, nVocabs] = size(barX);
assert(nDocs == nVocabs)
assert(max(barX(:)) <= 1)


for i=1:50
    [X, ~, ~] = prepareData('tdt2', 10, i, ...
        'dataType', 'tfidf', ...
        'modelType', 'P_CLEAN');
    fprintf('Index %d Data size %s \n', i, num2str(size(X)));
end
