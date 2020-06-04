function errors =  rmse_vs_speed_scale_factor(decoding_data, dt, n_bins_x, n_bins_y, n_electrodes, pixels_per_m1, tetrode_lamb_x, bin_size1, tetrode_z_coord, tetrode_lamb_ax, interspike_interval, sigmas, X_grid, Y_grid, x_real, y_real, step_count, n_tetrodes)


vel_const = [3:30]/10;
errors = zeros(1, size(vel_const, 2));
for i = 1:size(vel_const, 2);
    K = vel_const(i);
    disp(K);
    [log_post, decoded_position, n_enc_spikes, log_sum_lamb_ax_arr, dec_spikes] = compute_unsorted_posterior_position_precomp_v3(decoding_data, dt, n_bins_x, n_bins_y, n_tetrodes, n_electrodes, pixels_per_m1, tetrode_lamb_x, bin_size1, tetrode_z_coord, tetrode_lamb_ax, interspike_interval, sigmas, X_grid, Y_grid, step_count, K);
    rmse = xp_xr(decoded_position(:, 1), decoded_position(:, 2), x_real, y_real, 300);
    errors(i) = rmse;
    disp(errors);
    
end

plot(vel_const, errors);
xlabel('Proportionality Constant K');
ylabel('RMSE, cm');

end