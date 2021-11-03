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


    originalData = fea;
    if strcmp(Options.dataType, 'tfidf')
        fea = tfidf(fea, 0);
    elseif ~strcmp(Options.dataType, 'tf')
        error('NotImplement');
    end


    groundTruth = gnd;
    if strcmp(Options.modelType, 'C')
        X = fea;
        % originalData = fea;
        barX = X./sum(X, 1);
        
        IDidNormalized = true;
        assert(numel(unique(groundTruth)) == numTopics);
    elseif strcmp(Options.modelType, 'P_FOR_BASELINE')
        nDocs = size(fea, 1);
        X = (1/nDocs)*fea'*fea;
        % originalData = fea;
        IDidNormalized = true;
        barX = X;
    elseif strcmp(Options.modelType, 'P_PCA1')
        error('Deprecated!')

        if strcmp(source, 'tdt2')
            source = 'TDT2';
        elseif strcmp(source, 'reuters21578')
            source = 'reuters';
        else
            error ('mistypo');
        end

        filePCA = sprintf('%s/data/%s/%dClass/%d_pca.mat', filepath, ...
            source, numTopics, ind);
        if isfile(filePCA)
            fprintf('Load PCA input from mat \n');
            tmp = load(filePCA);
            X = tmp.X;
            barX = tmp.barX;
        else
            nDocs = size(fea, 1);
            P = (1/nDocs)*fea'*fea;

            X = P;
            barP = P ./ sum(P, 1);
            if isnan(sum(barP))
                error('Error in normalization');
            end

            [U, ~, ~] = svd(full(barP));
            barX = U(:, 1:2*numTopics)'*barP;
            save(filePCA, 'barX', 'X');
            fprintf('Save PCA to mat at %s \n', filePCA);
        end
        IDidNormalized = true;
        % originalData = fea;
    elseif strcmp(Options.modelType, 'P_PCA2')

        if strcmp(source, 'tdt2')
            source = 'TDT2';
        elseif strcmp(source, 'reuters21578')
            source = 'reuters';
        else
            error ('mistypo');
        end

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
            % [U, ~, ~] = svd(full(barP));
            % barX = U(:, 1:2*numTopics)'*barP;
            save(filePCA, 'barX', 'X');
            fprintf('Save PCA to mat at %s \n', filePCA);
        end

        IDidNormalized = true;

    elseif strcmp(Options.modelType, 'P_DIRTY')
        fea_t = fea;
        m = size(fea_t, 1);
        D_half = sum(fea_t*fea_t',2).^-0.5;
        doc_num = size(D_half, 1);
        D_half = spdiags(D_half, 0, doc_num, doc_num);
        fea_t = D_half * fea_t;
        P = (1/m)*fea_t'*fea_t;
        X = P;
    else
        error('NotImplement');
    end

    % zeroFilter = sum(X, 1) == 0;
    % if sum(zeroFilter) > 0
    %     if Options.verbose
    %         fprintf('WARNING: remove %d words \n', ...
    %             sum(full(zeroFilter)));
    %     end
    % end
    % X(:, zeroFilter) = [];
    % % TODO Im not sure about the next line
    % % originalData = originalData(:, zeroFilter);
    %
    % barX = X ./ sum(abs(X), 1);
    if ~IDidNormalized
        error('You forgot normlization')
    end
end

