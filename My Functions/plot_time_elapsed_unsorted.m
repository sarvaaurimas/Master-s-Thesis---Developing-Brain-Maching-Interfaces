function plot_time_elapsed_unsorted(timestamps, time_elapsed)
p = polyfit(timestamps, time_elapsed, 1);
best_fit = polyval(p,  timestamps);

plot(timestamps, time_elapsed, 'x');
hold on;
plot(timestamps, best_fit, '--');
legend('Real time values','Best fit line');
xlabel('N tetrodes used','FontSize',14);
ylabel('Decoding phase length, s','FontSize',14);
end