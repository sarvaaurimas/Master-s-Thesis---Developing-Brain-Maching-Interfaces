function plot_posterior_time_elapsed(cells, time_elapsed)
p = polyfit(cells, time_elapsed, 1);
best_fit = polyval(p,  cells);

plot(cells, time_elapsed, 'x');
hold on;
plot(cells, best_fit, '--');
legend('Real time values','Best fit line');
xlabel('Number of cells used','FontSize',14);
ylabel('Posterior generation time, s','FontSize',14);
end