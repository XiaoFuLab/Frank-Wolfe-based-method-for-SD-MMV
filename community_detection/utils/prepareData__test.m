%% Clear all things
clc; clear; close all; path(pathdef);

% [X, GT] = prepareData('DBLP', 2, 'trialInd', 0, 'isSecondOrder', 1);

[X, GT] = prepareData('YOUTUBE', 30, 'trialInd', 1);

[X, GT] = prepareData('EMAIL', 10, 'trialInd', 1);

[X, GT] = prepareData('DBLP', 1, 'isSecondOrder', false);

[X, GT] = prepareData('DBLP', 1);

size(X)
size(GT)

sumGT = sum(GT, 2);
norm(sum(GT, 2) - 1)
assert(norm(sum(GT, 2) - 1)<1e-4)

assert(sum(X - X', 'all') == 0)


numNodes = 0;
numComm = 0;
for i=1:5
    [X, GT] = prepareData('DBLP', i);
    numNodes = numNodes + size(X, 1);
    numComm = numComm + size(GT, 2);
end

[X, GT] = prepareData('DBLP', 8);
numNodesM1 = size(X, 1);

[X, GT] = prepareData('DBLP', 3);
numNodes1 = size(X, 1);

assert(numNodes - numNodes1 == numNodesM1)


