function rmse = plot_rmse(obsSpk, rms, dms, x_real, y_real, n, bin_size1, pixels_per_m1)
rmse = zeros([size(rms, 3)/10, 1]);
ncells = [1:size(rms, 3)/10]*10;
nbins = size(obsSpk, 2);
for k = 1:size(rms, 3)/10
    j = 10*k
    curr_rms = rms(:, :, 1:j);
    post = decode_calcBayesPost(obsSpk(1:j, :), curr_rms, dms, 1); % Compute posteriors for each step 
    max_post = reshape(max(max(post)), [nbins, 1]); % Calculate maximum posteriors in each step
    row = zeros(length(max_post), 1);
    col = zeros(length(max_post), 1);
    for i = 1:length(max_post)
        if isnan(max_post(i)) %for low cell sizes sometimes NaN values appear in post
            row(i) = row(i-1);
            col(i) = col(i-1);
        else
            [row(i), col(i)] = find(post(:, :, i) == max_post(i));% Get indexes of max posteriors in each step, sometimes throws errors if using too little cells
        end
    end
    y_pred = (row-0.5)*bin_size1/pixels_per_m1*100;
	x_pred = (col-0.5)*bin_size1/pixels_per_m1*100; % Predicted centers of bins in cms
    rmse(k, 1) = mean(sqrt((x_pred(1:n)-x_real(1:n)).^2 + (y_pred(1:n)-y_real(1:n)).^2));
end
plot(ncells, rmse);
title('RMSE vs number of cells');
xlabel('Number of Cells used');
ylabel('RMSE, cm');
