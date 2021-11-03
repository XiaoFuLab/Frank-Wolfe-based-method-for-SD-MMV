%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common/prox_ops')
addpath('~/code/matlab/common/PGD')

num_trials = 30;
list_hparams = 10:2:22;

success_rate = [];
for j=1:length(list_hparams)
    SNR = list_hparams(j);
    d = 50; L = 300; k = 40;

    success_count = 0;
    for i=1:num_trials
        fprintf('SNR=%d - Trial %d \n', SNR, i)
        A = rand(d, k);
        S = drchrnd(ones(1, k), L);
        S(1:k, 1:k) = eye(k);
        S(k+1:end, 1:k) = 0;
        X = A*S;
        noise = randn(size(X)); 
        sigma2 = sum(vecnorm(A*S, 2, 1).^2) / d / L / SNR;
        noise = sqrt(sigma2)*noise;
        X = X + noise;
        noise_eps = max(vecnorm(noise, 2, 1));
        noise_eps = noise_eps + 0; % (rand(1)/6.7 - .15)*noise_eps;
        delta = 4*noise_eps;

        lambda = greedy_pursuit(X, 'verbose', false, 'rank', k);
        if length(intersect(lambda, [1:k])) == k
            success_count = success_count + 1;
        end
    end
    success_rate(j) = success_count / num_trials;
end

figure();
plot(list_hparams, success_rate, '-o');
xlabel('SNR')
ylabel('Prob \Lambda = hat(\Lambda)')
% saveas(gcf, 'results/sd-somp-simulation.eps', 'epsc')

function r = drchrnd(a,n)
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
    r = r';
end

