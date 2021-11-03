%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')

addpath('./../../baselines/greedy_pursuit')
addpath('./../../baselines/FGNSR/matlab')
addpath('./../../../topic_modeling/baselines/PrecondSPA/')
addpath('./../../')
addpath('./../../topic_modeling/baselines/PrecondSPA')
addpath('./../../utils')
addpath('./../../fw_core')

%% UNIFORM A AND DIRICHLET S AND GAUSSIAN NOISE
M = 50; L = 200; N = 40;
num_trials = 100;
hparams = [5:0.5:16];
num_trials = 100;

success_count = [];

Tracking = struct;
Tracking.lambda = 1;
Tracking.M = M;
Tracking.L = L;
Tracking.N = N;
Tracking.num_trials = num_trials;
Tracking.hparams = hparams;

global_tic = tic;
output_dir = 'result/';
if ~isdir(output_dir)
    mkdir(output_dir)
    fprintf('Created output directory\n');
end

fprintf('Running X2\n')
num_exps = numel(hparams);
for exp_ind=1:num_exps
    fprintf('Exp %d/%d -  Hparam: %f\n', exp_ind, num_exps, ...
        hparams(exp_ind))
    Tracking.fw(exp_ind).success_count = [];
    Tracking.gp(exp_ind).success_count = [];
    Tracking.fg(exp_ind).success_count = [];
    Tracking.spa(exp_ind).success_count = [];

    for trial_ind=1:num_trials
        if mod(trial_ind, 10) == 0
            fprintf('Trial %d\n', trial_ind);
            fprintf('Total elapsed time: %.2f s\n', duration);
        end
            
        param = hparams(exp_ind);
        [X, pure_pixel_set] = generate_data(1, 3, M, N, L, 'SNR', param);

        % ==============================
        % GREEDY PUISUIT
        gp_tic = tic;
        gp_options = struct;
        gp_options.rank = N;
        gp_lambda_hat = greedy_pursuit(X, gp_options);
        Tracking.gp(exp_ind).success_count(trial_ind) = ...
            all(sort(gp_lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.gp(exp_ind).duration(trial_ind) = toc(gp_tic);

        % ==============================
        % FRANK-WOLFE
        fw_tic = tic;
        options = struct;
        options.maxIters = 100;
        options.verbose = 0;
        options.lambda = Tracking.lambda;
        options.debug = 0;
        options.N = N;
        options.backend = 'mex';
        options.dualGapThreshold=0;
        options.epsilon = 1e-5;
        [C_hat, fw_tracking] = fw(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), N, 1);

        Tracking.fw(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.fw(exp_ind).duration(trial_ind) = toc(fw_tic);
        Tracking.fw(exp_ind).lambda(trial_ind) = fw_tracking.lambda;


        % ==============================
        % FAST GRADIENT
        fg_tic = tic;
        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, N, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        Tracking.fg(exp_ind).success_count(trial_ind) = ...
            all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        Tracking.fg(exp_ind).duration(trial_ind) = toc(fg_tic);


        % ==============================
        % SPA
        i_tic = tic;
        K = SPAselect(X, N);
        Tracking.spa(exp_ind).success_count(trial_ind) = ...
            all(sort(K(:)) == sort(pure_pixel_set(:)));
        Tracking.spa(exp_ind).duration(trial_ind) = toc(i_tic);

        duration = toc(global_tic);
        Tracking.duration = duration;

        % ==============================
        % FRANK-WOLFE
        fw_tic = tic;
        options = struct;
        options.maxIters = 100;
        options.verbose = 0;
        options.lambda = 0;
        options.debug = 0;
        options.N = N;
        options.backend = 'mex';
        options.dualGapThreshold=0;
        options.epsilon = 1e-5;
        [C_hat, fw_tracking] = fw(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), N, 1);

        Tracking.fw_(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.fw_(exp_ind).duration(trial_ind) = toc(fw_tic);
        Tracking.fw_(exp_ind).lambda(trial_ind) = fw_tracking.lambda;


        save(sprintf('%s/result.mat', output_dir), 'Tracking')
    end
end

fw = arrayfun(@(x) mean(x.success_count), Tracking.fw);
gp = arrayfun(@(x) mean(x.success_count), Tracking.gp);
fg = arrayfun(@(x) mean(x.success_count), Tracking.fg);
spa = arrayfun(@(x) mean(x.success_count), Tracking.spa);
