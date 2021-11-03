%% Clear all things
clc; clear; close all; path(pathdef);

for i=1:5
    [XO, GTO] = prepareData('MAG', i);   
    L = size(XO, 1);
    for j=1:50
        sampleIndex = randsample(1:L, 5000);
        X = XO(sampleIndex, sampleIndex);
        GT = GTO(sampleIndex, :);
        save(sprintf('./../data/coauthorship/random-split2/MAG%d_trial_%d.mat', i, j), ...
            'X', 'GT');
    end
    i
end

