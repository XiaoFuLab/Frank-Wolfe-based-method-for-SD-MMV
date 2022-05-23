function [barX, groundTruth, X, originalData] = prepareData(source, numTopics, ind,...
        varargin) 
% -----------------------------------------------------------------------
% 
% X: V x num_docs
% ground_truth: num_docs. Values are topic indices of documents X


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 1);
    p.addOptional('debug', 0);
    p.addOptional('dataType', '');
    p.addOptional('modelType', '');
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    Options = p.Results;

    IDidNormalized = false;
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    if strcmp(source, 'reuters21578')
        Data = load([filepath, '/data/Reuters21578.mat']);
        DataAux = load([filepath, ...
            sprintf('/data/%dClass/%d.mat', numTopics, ind)]);
        fea = Data.fea(DataAux.sampleIdx,:);
        fea(:, DataAux.zeroIdx) = [];
        gnd = Data.gnd(DataAux.sampleIdx,:);
    elseif strcmp(source, 'tdt2')
        Data = load([filepath, '/data/TDT2/TDT2.mat']);
        DataAux = load([filepath, ...
            sprintf('/data/TDT2/%dClass/%d.mat', numTopics, ind)]);
        fea = Data.fea(DataAux.sampleIdx,:);
        fea(:, DataAux.zeroIdx) = [];
        gnd = Data.gnd(DataAux.sampleIdx,:);
    elseif strcmp(source, '20newsgroups')
        Data = load([filepath, '/data/20newsgroups/20Newsgroups.mat']);
        fea = Data.fea;
        gnd = Data.gnd;

        assert(max(gnd) == numel(unique(gnd)));
        topicFilter = ismember(gnd, [1:numTopics]);
        fea = fea(topicFilter, :);
        gnd = gnd(topicFilter);
    else
        error('NotImplement')
    end


    fea = tfidf(fea, 0);


    groundTruth = gnd;

    if strcmp(source, 'tdt2')
        source = 'TDT2';
    elseif strcmp(source, 'reuters21578')
        source = 'reuters';
    else
        error ('mistypo');
    end

    originalData = fea;
    filePCA = sprintf('%s/data/%s/%dClass/%d_pca.mat', filepath, ...
        source, numTopics, ind);

    zeroFilter = sum(fea) == 0;
    originalData(:, zeroFilter) = [];

    if isfile(filePCA)
        fprintf('Load PCA input from mat \n');
        tmp = load(filePCA);
        X = tmp.X;
        barX = tmp.barX;
    else
        nDocs = size(fea, 1);
        fea(:, zeroFilter) = [];
        P = (1/nDocs)*fea'*fea;

        X = P;
        barP = P ./ sum(P, 1);
        if isnan(sum(barP(:)))
            error('Error in normalization');
        end

        [U, ~, ~] = svds(barP, 2*numTopics);
        barX = U(:, 1:2*numTopics)'*barP;
        save(filePCA, 'barX', 'X');
        fprintf('Save PCA to mat at %s \n', filePCA);
    end

    IDidNormalized = true;

    if ~IDidNormalized
        error('You forgot normlization')
    end
end

