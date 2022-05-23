% This script produces results in Fig 4a in section V.A.
% You might want to change some settings at lines 15-16.


% ---------------------------------------------------------
%% Clear all things
clc; clear; close all; path(pathdef);
addpath('../baselines/PrecondSPA/')
addpath('../baselines/FGNSR/matlab/')
addpath('../utils')
addpath('../fw_core')
addpath('../')

% This is a courser result. 
hparams = [4:1:16];
num_trials = 5;
% For exact results shown in Fig 4a, please use 
% hparams = [4:1:16];
% num_trials = 50;


M=50; N=200; K=40;
alpha = ones(K, 1);
Tracking = struct;
Tracking.M = M;
Tracking.N = N;
Tracking.K = K;
Tracking.num_trials = num_trials;
Tracking.hparams = hparams;


fprintf('Running X2\n')
num_exps = numel(hparams);
for exp_ind=1:num_exps
    fprintf('Exp %d/%d -  Hparam: %f\n', exp_ind, num_exps, ...
        hparams(exp_ind))
    Tracking.fw(exp_ind).success_count = [];
    Tracking.fg(exp_ind).success_count = [];
    Tracking.spa(exp_ind).success_count = [];

    for trial_ind=1:num_trials
        if mod(trial_ind, 10) == 0
            fprintf('Trial %d\n', trial_ind);
            fprintf('Total elapsed time: %.2f s\n', duration);
        end
            
        SNR = hparams(exp_ind);

        W = rand(M, K);
        H = zeros(K, N);
        H(:, 1:K) = eye(K);
        H(:, K+1:end) = dirichlet_rnd(alpha, N-K);
        Y = W*H;
        SNR = 10^(SNR/10);
        noise = randn(size(Y)); 
        sigma2 = sum(vecnorm(Y, 2, 1).^2) / M / N / SNR;
        noise = sqrt(sigma2)*noise;
        X = Y + noise;
        indices = randperm(N);
        X = X(:, indices);
        H = H(:, indices);
        r_pure_pixel_set = [];
        pure_pixel_set = 1:K;
        for ii=1:numel(pure_pixel_set)
            r_pure_pixel_set(end+1) = find(indices == pure_pixel_set(ii));
        end
        pure_pixel_set = r_pure_pixel_set;


        % ==============================
        % FAST GRADIENT
        fg_tic = tic;
        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        Tracking.fg(exp_ind).success_count(trial_ind) = ...
            all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        Tracking.fg(exp_ind).duration(trial_ind) = toc(fg_tic);

        % ==============================
        % SPA
        i_tic = tic;
        Lambda_hat_SPA = SPAselect(X, K);
        Tracking.spa(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat_SPA(:)) == sort(pure_pixel_set(:)));
        Tracking.spa(exp_ind).duration(trial_ind) = toc(i_tic);

        % ==============================
        % FRANK-WOLFE
        fw_tic = tic;
        options = struct;
        options.epsilon = 1e-5;
        W_hat_spa = X(:, Lambda_hat_SPA);
        H_hat_spa = constrained_ls(W_hat_spa, X, 'max_iters', 100);
        %   Use SPA solution to determine lambda
        options.lambda = norm(X - W_hat_spa*H_hat_spa, 'fro')/K;

        options.maxIters = 100;
        options.verbose = 0;
        options.debug = 0;
        options.backend = 'mex';
        [C_hat, fw_tracking] = MERIT(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);

        Tracking.fw(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.fw(exp_ind).duration(trial_ind) = toc(fw_tic);
        Tracking.fw(exp_ind).lambda(trial_ind) = fw_tracking.lambda;
    end
end

fprintf('Running X2\n')
fw = arrayfun(@(x) mean(x.success_count), Tracking.fw);
fg = arrayfun(@(x) mean(x.success_count), Tracking.fg);
spa = arrayfun(@(x) mean(x.success_count), Tracking.spa);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineMarkerSize', 9)
set(groot, 'DefaultLineLineWidth', 1.4)
set(groot,'defaultAxesFontSize',20)

fw_color = [0 0.4470 0.7410];
fg_color = '#EDB120';
spa_color = [0.4940, 0.1840, 0.5560];

figure('DefaultAxesFontSize', 18);
plot(hparams, fw, '-s',  'DisplayName', '\texttt{MERIT}', 'Color', fw_color)
hold on
plot(hparams, fg, '-o',  'DisplayName', '\texttt{FastGradient}', 'Color', fg_color)
plot(hparams, spa, '-x',  'DisplayName', '\texttt{SPA}', 'Color', spa_color)

xlabel('SNR');
ylabel('success rate');
legend('Location', 'northwest', 'Interpreter', 'latex');
axis tight;

