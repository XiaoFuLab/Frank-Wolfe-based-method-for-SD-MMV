% Produce experiment results in Fig 2, and Fig 3 in Section V.A.
% You might want to change a setting in line 13.

% --------------------------------------------------------------
%% Clear all things
clc; clear; close all; path(pathdef);
addpath('../baselines/PrecondSPA/')
addpath('../baselines/FGNSR/matlab/')
addpath('../utils')
addpath('../fw_core')
addpath('../')

M=50; N=55; K=10;
% This is a coarser version. For exact result as shown in the paper, use list_SNR = [5:1:30];
list_SNR = [5:3:30];

RAE_fw = [];
RAE_spa = [];
RAE_fg = [];

succ_fw = [];
succ_spa = [];
succ_fg = [];

MRSA_fw = [];
MRSA_spa = [];
MRSA_fg = [];

t0 = tic;
for exp_ind=1:numel(list_SNR)
    fprintf('%d/%d - Elapsed time: %f \n', exp_ind, numel(list_SNR), toc(t0));

    for trial_ind=1:20
        W = rand(M, K);
        W = W./sum(W, 1);

        H = zeros(K, N);
        H(:, 1:K) = eye(K);
        col = K+1;
        for i=1:K
            for ii=i+1:K
                H(i, col) = 0.5;
                H(ii, col) = 0.5;
                col = col+1;
            end
        end
        Y = W*H;
        pure_pixel_set = 1:K;

        w_bar = mean(W, 2);
        original_noise = zeros(M, N);
        original_noise(:, K+1:end) = Y(:, K+1:end) - w_bar;

        noise_energy = sqrt(norm(Y, 'fro')^2/10^(list_SNR(exp_ind)/10));
        noise = original_noise/norm(original_noise, 'fro')*noise_energy;
        X = Y + noise;
        SNR(exp_ind, trial_ind) = 10*log10(norm(Y, 'fro')^2/norm(noise, 'fro')^2);

        new_order = randperm(N);
        r_pure_pixel_set = [];
        for i=1:numel(pure_pixel_set)
            r_pure_pixel_set(end+1) = find(new_order == pure_pixel_set(i));
        end
        pure_pixel_set = r_pure_pixel_set;
        X = X(:, new_order);

        % ------------- SPA ------------- 
        Lambda_hat = SPAselect(X, K);
        succ_spa(exp_ind, trial_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        W_hat_spa = X(:, Lambda_hat);
        H_hat_spa = constrained_ls(W_hat_spa, X, 'max_iters', 100);
        RAE_spa(exp_ind, trial_ind) = 1-norm(X - W_hat_spa*H_hat_spa, 'fro')/norm(X, 'fro');
        MRSA_spa(exp_ind, trial_ind) = compute_MRSA(W_hat_spa, W);


        % ------------- FRANK-WOLFE ------------- 
        fw_tic = tic;
        options = struct;
        %  Use SPA solution to estimate lambda
        options.lambda = norm(X - W_hat_spa*H_hat_spa, 'fro')/K;
        options.maxIters = 100;
        options.verbose = 0;
        options.debug = 0;
        options.N = K;
        options.backend = 'mex';
        options.dualGapThreshold=0;
        options.epsilon = 1e-5;
        [C_hat, fw_tracking] = MERIT(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);

        succ_fw(exp_ind, trial_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));

        W_hat = X(:, Lambda_hat);
        H_hat = constrained_ls(W_hat, X, 'max_iters', 100);
        RAE_fw(exp_ind, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        MRSA_fw(exp_ind, trial_ind) = compute_MRSA(W_hat, W);


        % ------------- FastGradient ------------- 
        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        succ_fg(exp_ind, trial_ind) = all(sort(K_fg(:)) == sort(pure_pixel_set(:)));

        W_hat = X(:, K_fg);
        H_hat = constrained_ls(W_hat, X, 'max_iters', 100);
        RAE_fg(exp_ind, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        MRSA_fg(exp_ind, trial_ind) = compute_MRSA(W_hat, W);
    end
end


SNR = list_SNR;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineMarkerSize', 9)
set(groot, 'DefaultLineLineWidth', 1.4)
set(groot, 'DefaultAxesFontSize', 14);

figure();
plot(SNR, mean(RAE_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}')
hold on
plot(SNR, mean(RAE_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}')
plot(SNR, mean(RAE_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}')
xlabel('SNR')
ylabel('\texttt{Relative approximation error}', 'Interpreter',  'latex')
legend('Location',  'NorthWest')

tmp = max([max(MRSA_spa(:)) max(MRSA_fw(:)) max(MRSA_fg(:))]);
figure();
plot(SNR, mean(MRSA_fw, 2)/tmp*100, '-s', 'DisplayName',  '\texttt{MERIT}')
hold on
plot(SNR, mean(MRSA_fg, 2)/tmp*100, '-o', 'DisplayName',  '\texttt{FastGradient}')
plot(SNR, mean(MRSA_spa, 2)/tmp*100, '-x', 'DisplayName',  '\texttt{SPA}')
ylabel('\texttt{MRSA}', 'Interpreter',  'latex')
xlabel('SNR')
legend

figure();
plot(SNR, mean(succ_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}')
hold on
plot(SNR, mean(succ_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}')
plot(SNR, mean(succ_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}')
ylabel('\texttt{success rate}', 'Interpreter',  'latex')
xlabel('SNR')
legend('Location',  'NorthEast')

