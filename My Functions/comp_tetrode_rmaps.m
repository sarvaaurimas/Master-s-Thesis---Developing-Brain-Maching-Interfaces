function tetrode_rate_maps = comp_tetrode_rmaps(encoding_vectors, n_tetrode_spikes, mu, pixels_per_m1, max_time, n_bins_y, n_bins_x, xy1, X_bin, Y_bin, B_x)
tetrode_rate_maps = zeros(n_bins_y, n_bins_x, size(n_tetrode_spikes, 1));
for i = 1:(size(n_tetrode_spikes)-1)
    disp(i);
    if i == 1
        [lamb_x1, pi_x] = lambda_x(X_bin, Y_bin, B_x, encoding_vectors(1:n_tetrode_spikes(1), 2:3), mu, pixels_per_m1, max_time/4, n_bins_y, n_bins_x, xy1);
    else
        [lamb_x1, pi_x] = lambda_x(X_bin, Y_bin, B_x, encoding_vectors(sum(n_tetrode_spikes(1:i)):sum(n_tetrode_spikes(1:i+1)), 2:3), mu, pixels_per_m1, max_time/4, n_bins_y, n_bins_x, xy1);
    end
    tetrode_rate_maps(:, :, i) = lamb_x1;
    surf(X_bin, Y_bin, lamb_x1);
    pause;
end
end