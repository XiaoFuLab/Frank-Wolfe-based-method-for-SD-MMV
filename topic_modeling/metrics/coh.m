%% Clear all things
clc; clear; close all; path(pathdef);

load('./../data/Reuters21578.mat')
fea;
V = size(fea, 2);

topic_a = [1 10 20];
topic_docs = fea(topic_a, :);

num_docs = size(topic_docs, 1);
topic_docs = mat2cell(topic_docs, ones(1, num_docs), [V]);
score = 0;
for i=1:V
    for j=1:V
        freq_v2 = sum(cellfun(@(x) ismember(j, x), topic_docs), 'all');
        freq_v1v2 = sum(cellfun(@(x) ismember(j, x) & ismember(i, x),...
            topic_docs), 'all');
        score = score + log((freq_v1v2+0.01) / freq_v2);
    end
end

