function [rmse, err] = xp_xr(x_pred, y_pred, x_real, y_real, n)
%[rmse, err] = xp_xr(x_pred, y_pred, x_real, y_real, n, max_post)
err = sqrt((x_pred(1:n)-x_real(1:n)).^2 + (y_pred(1:n)-y_real(1:n)).^2);
rmse = mean(err);

% % Top comparison plot
% subplot(3, 1, 1);
% plot(x_pred(1:n)); %Function for plotting x vs xpred for 
% hold on;
% plot(x_real(1:n, 1)); 
% legend('X_{predicted}','X_{real}'); 
% xlabel('Time, s');
% ylabel('X, cm');
% 
% subplot(3, 1, 2);
% plot(y_pred(1:n)); %Function for plotting x vs xpred for 
% hold on;
% plot(y_real(1:n, 1)); 
% legend('Y_{predicted}','Y_{real}'); 
% xlabel('Time, s');
% ylabel('Y, cm');
% 
% 
% %Bottom rmse plot
% subplot(3, 1, 3);
% plot(err(1:n));
% %yyaxis right
% %plot(max_post(1:n));
% xlabel('Time, s');
% ylabel('RMSE, cm');
% title(['Mean RMSE = ', num2str(rmse), 'cm']);

