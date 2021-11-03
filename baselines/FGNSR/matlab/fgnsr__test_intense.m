%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common/FW')
addpath('~/code/matlab/common/prox_ops')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/')

addpath('./../../')

d = 1000; L = 20000; k = 30;
A = rand(d, k);
S = dirichlet_rnd(ones(1, k), L);

pure_pixel_set = randsample(L, k)';
for i=1:k
    S(:, pure_pixel_set(i)) = zeros(k, 1);
    S(i, pure_pixel_set(i)) = 1;
end

X = A*S;
SNR = 10; SNR = 10^(SNR/10);
noise = randn(size(X)); 
sigma2 = sum(vecnorm(A*S, 2, 1).^2) / d / L / SNR;
noise = sqrt(sigma2)*noise;
% X = X + noise;
tic
[C, K, ~, ~, ~, ~, t] = fgnsr(X, k, 'maxiter', 10, 'debug', 0);
toc

% disp('True basic indices: ')
% disp(sort(pure_pixel_set(:))');
% disp('Computed indices:')
% disp(sort(K(:))');
%
% disp('==========================');
% disp('done.');
% 
% figure();
% semilogy(t.time, t.rel_approx_error);
% xlabel('time (second)')
% ylabel('Relative approx error')
% 
% figure();
% semilogy(t.rel_approx_error);
% xlabel('Iter')
% ylabel('Relative approx error')
% 
% figure();
% plot(t.bytes/1024/1024);
% xlabel('iter');
% ylabel('Mb');
% title(sprintf('224x%d (10)', L))
% saveas(gcf, sprintf('./../../results/m-fg-%d.eps', L), 'epsc')
