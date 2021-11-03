%% Clear all things
clc; clear; close all; path(pathdef);

[X, GT] = prepareData('DBLP', 5, 'trialInd', 0);
N = size(GT, 2);
[V, D] = eigs(X, N);
V = V';

figure()
scatter(1:size(V, 2), vecnorm(V, 2, 1));

colEnergy = sort(vecnorm(V, 1, 1), 'descend');
totalEnergy = sum(colEnergy);
colEnergy = colEnergy;

threshold = 0.7*totalEnergy;
L = size(X, 2);
currentSum = 0;
for i=1:L
    currentSum = currentSum + colEnergy(i);
    if currentSum > threshold
        colEnergy(i)
        break
    end
end
