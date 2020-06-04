function plot_rmse_vs_n_vel_steps(errors, steps)
plot(steps, errors);
xlabel('N previous bins used');
ylabel('RMSE, cm');
end