%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')

addpath('./../../baselines/greedy_pursuit')
addpath('./../../baselines/FGNSR/matlab')
addpath('./../../')
addpath('./../../topic_modeling/baselines/PrecondSPA')
addpath('./../../utils')
addpath('./../../fw_core')

%% UNIFORM A AND DIRICHLET S AND GAUSSIAN NOISE
global_tic = tic;
output_dir = './result/';
load(sprintf('%s/result.mat', output_dir))

fprintf('Running X2\n')
fw = arrayfun(@(x) mean(x.success_count), Tracking.fw);
fg = arrayfun(@(x) mean(x.success_count), Tracking.fg);
spa = arrayfun(@(x) mean(x.success_count), Tracking.spa);
fw_ = arrayfun(@(x) mean(x.success_count), Tracking.fw_);

plot_success_rate(Tracking.hparams, {fw, fw_, fg, spa, }, ...
    {'FW', 'FW(0)', '\texttt{FastGradient}', '\texttt{SPA}', }, ...
    output_dir, 'SNR', 'export', 1, 'ylabel', 'success rate');


% fw = arrayfun(@(x) mean(x.lambda), Tracking.fw);
% figure();
% plot(Tracking.hparams, fw, 'DisplayName', 'fw')
% hold on
% legend
