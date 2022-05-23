% Display the results on a plot

colors{1} = 'b-+' ; % SPA
colors{2} = 'r' ; % SDP-SPA
colors{3} = 'co' ; % PW-SPA
colors{4} = 'k--' ; % SPA-SPA
colors{5} = 'g-x' ; % VCA
colors{6} = 'm-d' ; % XRAY

figure; 
for i = 1 : nalgo 
    plot(error, results(i,:), colors{i}, 'LineWidth', 2, 'MarkerSize',8); hold on;
end
set(gca,'FontSize',14); 
xlabel('Noise level (\delta)','FontSize',16); ylabel('Percentage of columns of W exctracted','FontSize',16) 
legend('SPA','SDP-SPA','PW-SPA','SPA-SPA','VCA','XRAY'); 