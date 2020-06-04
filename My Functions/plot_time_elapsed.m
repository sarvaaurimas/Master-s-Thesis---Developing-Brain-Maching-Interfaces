function plot_time_elapsed(timestamps, time_elapsed)
timestamps = timestamps / 50 / 60; % 50 is the sampling freq of timestamps and 60 to convert s -> min
p = polyfit(timestamps, time_elapsed, 1);
best_fit = polyval(p,  timestamps);

plot(timestamps, time_elapsed, 'x');
hold on;
plot(timestamps, best_fit, '--');
legend('Real time values','Best fit line');
xlabel('Experiment length, min','FontSize',14);
ylabel('Ratemap generation time, s','FontSize',14);
end