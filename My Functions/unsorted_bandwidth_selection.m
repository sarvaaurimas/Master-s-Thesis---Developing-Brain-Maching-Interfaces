function [bw_opt, rmse_grid, h_a_grid, h_x_grid] = unsorted_bandwidth_selection(encoding_vectors, n_tetrode_spikes, sigmas, encoding_time, xy1, dt, dec_raw_data, n_bins_x, n_bins_y, n_bins_z, n_tetrodes, n_electrodes, interspike_interval, X_bin, Y_bin, pixels_per_m1, bin_size1, x_real, y_real)

h_as = [5, 10, 15, 20, 25, 30];
h_xs = [2.5, 5, 7.5, 10];
h_xs = [2.5, 5, 7.5];
h_as = [7.5];
[h_x_grid, h_a_grid] = meshgrid(h_xs, h_as);

rmse_grid = zeros(size(h_as, 2), size(h_xs, 2));
disp(size(rmse_grid));
for i = 1:size(h_as, 2)
    h_a = h_as(i);
    for j = 1:size(h_xs, 2)
        disp([i, j]);
        h_x = h_xs(j);
        B_x = [h_x, h_x];
        B_dec = [h_a, h_x, h_x];
        disp(B_dec);
        pi_x = occupancy_x(X_bin, Y_bin, B_x, pixels_per_m1, encoding_time, n_bins_y, n_bins_x, xy1);
        tetrode_mu = n_tetrode_spikes/encoding_time; % Array of inidividual tetrode mu's
        tetrode_lamb_x = tetrode_lambda_x(X_bin, Y_bin, B_x, n_bins_y, n_bins_x, encoding_vectors, n_tetrode_spikes, pi_x, n_tetrodes, tetrode_mu);
        [tetrode_lamb_ax, X_bin_3d, Y_bin_3d, Z_bin_3d, tetrode_z_coord] = compute_tetrode_lambda_ax_v2(B_dec, encoding_vectors, tetrode_mu, pixels_per_m1, n_bins_y, n_bins_x, n_bins_z, pi_x, bin_size1, n_tetrodes, n_tetrode_spikes);

        [log_post, decoded_position, n_dec_spikes, log_sum_lamb_ax_arr, dec_spikes] = compute_unsorted_posterior_position_precomp_v2(...
                                                           dec_raw_data, dt, n_bins_x, n_bins_y, n_tetrodes, ...
                                                           n_electrodes, pixels_per_m1, tetrode_lamb_x, bin_size1, ...
                                                           tetrode_z_coord, tetrode_lamb_ax, ...
                                                           interspike_interval, sigmas);

        [rmse, err] = xp_xr(decoded_position(:, 1), decoded_position(:, 2), x_real, y_real, size(decoded_position, 1));
        rmse_grid(i, j) = rmse;
        disp(rmse);
    end
end

min_rmse = min(min(rmse_grid));
[ha_indeks, hx_indeks] = find(rmse_grid == min_rmse);
ha_opt = h_as(ha_indeks);
hx_opt = h_xs(hx_indeks);
bw_opt = [ha_opt, hx_opt, hx_opt];

end