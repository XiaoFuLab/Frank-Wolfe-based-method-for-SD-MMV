%% Clear all things
clc; clear; close all; path(pathdef);

FW = [ 775272 774276 786028 789212 787968 797856 808376 800608 809228 808716] - 740000;
FW_ = [774688 784184 782600 785176 791368 791560 807376 816212 805776 803452] - 740000;

FG = [ 848980 1057124 1317024 1711348 2217360 2830624 3553624 4378504 5319204 6369240] - 740000;
FW = FW/1024/1024;
FW_ = FW_/1024/1024;
FG = FG/1024/1024;

L = 1000:1000:10000;


% ==============================
% BEGIN: RRS
% ------------------------------

figure();
linewidth = 1.4;
markersize = 9;
marker_pools = {'s', 'o', '^', 'x', 's', 'd', '^', 'v'};
linestyle_pools = {'-', '--', '-.'};

lines{1} = FW;
lines{2} = FW_;
lines{3} = FG;
names = {'FW', 'FW(0)', '\texttt{FastGradient}', };
for i=1:numel(lines)
    semilogy(L, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
        'DisplayName', names{i}, ...
        'LineWidth', linewidth,  ...
        'MarkerSize', markersize ...
        );
    hold on
end
xlabel('N');
yticks([0.1 0.30 1.0 5.0 10.0])
ylabel('Memory cost in RSS (GB)');
legend('Location', 'northwest', 'Interpreter', 'latex');
axis tight;
set(gca, 'FontSize', 15);
output_dir = './results/';
path_to_file = sprintf('%s/mem.eps', output_dir);

% exportgraphics(gcf, path_to_file, 'resolution', 300);
saveas(gcf, path_to_file, 'epsc')

fprintf(sprintf('Exported to %s\n', path_to_file));
% ------------------------------
% END: RRS
% ==============================



% ==============================
% BEGIN: NNZ
% ------------------------------
% FW = [];
% FG = [];
% for i=1:numel(L)
%     t = load(sprintf('./results/memory_%d/FW/result.mat', L(i)));
%     FW(i) = nnz(t.Tracking.trackingMethod{1}.hatC);
%
%     t = load(sprintf('./results/memory_%d/FG/result.mat', L(i)));
%     FG(i) = t.Tracking.trackingMethod{1}.nnz(end);
% end
%
% figure();
%
%
% linewidth = 1.4;
% markersize = 9;
% marker_pools = {'s', 'o', '^', 'x', 's', 'd', '^', 'v'};
% linestyle_pools = {'-', '--', '-.'};
%
%
% lines{1} = FW;
% lines{2} = FG;
% names = {'FW', 'FastGradient'};
% for i=1:numel(lines)
%     semilogy(L, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
%         'DisplayName', names{i}, ...
%         'LineWidth', linewidth,  ...
%         'MarkerSize', markersize ...
%         );
%     hold on
% end
% xlabel('N');
% ylabel('nnz');
% legend('Location', 'southeast');
% %
% set(gca, 'FontSize', 15);
% output_dir = './results/';
% path_to_file = sprintf('%s/nnz.eps', output_dir);
%
% % exportgraphics(gcf, path_to_file, 'resolution', 300);
% saveas(gcf, path_to_file, 'epsc')
%
% fprintf(sprintf('Exported to %s\n', path_to_file));

% ------------------------------
% END: NNZ
% ==============================

