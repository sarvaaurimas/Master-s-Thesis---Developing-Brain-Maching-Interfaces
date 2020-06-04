function [rmse, err] = pred_real_trajectories(x_pred, y_pred, x_real, y_real, n)
%[rmse, err] = xp_xr(x_pred, y_pred, x_real, y_real, n, max_post)
err = sqrt((x_pred(1:n)-x_real(1:n)).^2 + (y_pred(1:n)-y_real(1:n)).^2);
rmse = mean(err);

plot(x_pred(1:n), y_pred(1:n)); %Function for plotting x vs xpred for 
hold on;
plot(x_real(1:n, 1), y_real(1:n, 1)); 
legend('X_{predicted}','X_{real}'); 
xlabel('x, cm');
ylabel('y, cm');
title(['Mean RMSE = ', num2str(rmse), 'cm']);

