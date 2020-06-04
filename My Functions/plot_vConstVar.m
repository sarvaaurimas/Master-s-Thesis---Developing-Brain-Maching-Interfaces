function rmse = plot_vConstVar(n_test_bins, obsSpk, rms, dms, x_real, y_real, bin_size1, pixels_per_m1, tBin)
size = 60;
rmse = zeros([1, size]);
vConstArr = (1:size)/20; % increments of 0.167
for j = (1:size)
    j
    vConst = vConstArr(j);
    post = decode_calcBayesPost(obsSpk(:, 1:n_test_bins), rms, dms, tBin, vConst); % Compute posteriors for each step 
    max_post = reshape(max(max(post)), [n_test_bins, 1]); % Calculate maximum posteriors in each step
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
    rmse(1, j) = mean(sqrt((x_pred(1:n_test_bins)-x_real(1:n_test_bins)).^2 + (y_pred(1:n_test_bins)-y_real(1:n_test_bins)).^2));
end

plot(vConstArr, rmse);
title('RMSE vs K');
xlabel('K');
ylabel('RMSE, cm');