function [out] = plot_success_rate(hparams, lines, names, output_dir, xlabel_name, varargin) 
% -------------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -------------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('fig_visible', 1);
    p.addOptional('export', 1);
    p.addOptional('location', 'southeast');
    p.addOptional('title', '');
    p.addOptional('marker_pools', {'s', 'o', '^', 'x', 's', 'd', '^', 'v'});
    p.addOptional('linestyle_pools', {'-', '--', '-.', '-'});
    p.addOptional('ylabel', '');

    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    figure('visible', options.fig_visible);
    linewidth = 1.4;
    markersize = 9;
    marker_pools = options.marker_pools;
    linestyle_pools = options.linestyle_pools;
    assert(numel(marker_pools) >= numel(lines))
    for i=1:numel(lines)
        plot(hparams, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
            'DisplayName', names{i}, ...
            'LineWidth', linewidth,  ...
            'MarkerSize', markersize ...
            );

        hold on
    end
    xlabel(xlabel_name);
    ylabel(options.ylabel);
    legend('Location', options.location, 'Interpreter', 'latex');

    title(options.title);
    % ylim([-0.1 1.1]);
    axis tight;
    set(gca, 'FontSize', 15);
    if options.export
        path_to_file = sprintf('%s/succ_rate.eps', output_dir);

        % exportgraphics(gcf, path_to_file, 'resolution', 300);
        saveas(gcf, path_to_file, 'epsc')

        fprintf(sprintf('Exported to %s\n', path_to_file));
    end
end
