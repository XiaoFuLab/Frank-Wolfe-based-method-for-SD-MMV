%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common/prox_ops')
addpath('~/code/matlab/common/PGD')

num_trials = 30;
hparams = [3];
success_rate = [];
for exp_ind=1:length(hparams)
    SNR = hparams(exp_ind);
    M = 50; L = 55; N = 10;

    success_count = 0;
    for trial_ind=1:num_trials
        fprintf('epsilon=%f - Trial %d \n', hparams(exp_ind), trial_ind)
        A = rand(M, N);
        S = zeros(N, L);
        S(:, 1:N) = eye(N);
        t=N+1;
        for i=1:N
            for j=i+1:N
                S(:, t) = 0.5*(S(:, i) + S(:, j));
                t = t+1;
            end
        end
        pure_pixel_set = [1:N];
        X = A*S;
        Noise = zeros(size(X));
        Noise(:, 1:N) = 0;
        center = mean(A, 2);
        Noise(:, N+1:end) = X(:, N+1:end)-center;
        epsilon = hparams(exp_ind);
        Noise = epsilon*Noise/norm(Noise, 'fro');
        X = X + Noise;
        
        lambda = greedy_pursuit(X, 'verbose', false, 'rank', N);
        if length(intersect(lambda, [1:N])) == N
            success_count = success_count + 1;
        end
    end
    success_rate(exp_ind) = success_count / num_trials;
end

figure();
plot(hparams, success_rate, '-o');
xlabel('epsilon')
ylabel('Prob \Lambda = hat(\Lambda)')
% saveas(gcf, 'results/sd-somp-simulation.eps', 'epsc')

function r = drchrnd(a,n)
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
    r = r';
end

