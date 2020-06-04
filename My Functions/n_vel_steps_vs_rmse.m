function errors = n_vel_steps_vs_rmse(decoding_data, dt, n_bins_x, n_bins_y, n_tetrodes, n_electrodes, pixels_per_m1, tetrode_lamb_x, bin_size1, tetrode_z_coord, tetrode_lamb_ax, interspike_interval, sigmas, X_grid, Y_grid, x_real, y_real)

n_steps = [20:30];
errors = zeros(1, size(n_steps, 2));
for i = 1:size(n_steps, 2);
    step_count = n_steps(i);
    disp(step_count);
    [log_post, decoded_position, n_enc_spikes, log_sum_lamb_ax_arr, dec_spikes] = compute_unsorted_posterior_position_precomp_v3(decoding_data, dt, n_bins_x, n_bins_y, n_tetrodes, n_electrodes, pixels_per_m1, tetrode_lamb_x, bin_size1, tetrode_z_coord, tetrode_lamb_ax, interspike_interval, sigmas, X_grid, Y_grid, step_count);
    rmse = xp_xr(decoded_position(:, 1), decoded_position(:, 2), x_real, y_real, 300);
    errors(i) = rmse;
    disp(errors);
    
end